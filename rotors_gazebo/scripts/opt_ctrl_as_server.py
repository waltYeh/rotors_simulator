#!/usr/bin/env python
from casadi import *
import casadi.tools as ctools
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from trajectory_msgs.msg import MultiDOFJointTrajectory, MultiDOFJointTrajectoryPoint
from geometry_msgs.msg import Vector3
from geometry_msgs.msg import Twist
from geometry_msgs.msg import Transform, Quaternion
import rospy
from rotors_gazebo.srv import *
import tf
import sys
def project_parameters():
	kF=10
	kM=2
	mass=1.0
	Len=0.225
	gravity=9.81
	return gravity, mass, kF, kM, Len
def diff_eq_time_trans(tau,x,u,t_trans):
	gravity, mass, kF, kM, Len = project_parameters()
	dangvel_Body_dt = vertcat(kF*Len*(u[1]-u[3]), kF*Len*(-u[0]+u[2]), kM*(u[0]-u[1]+u[2]-u[3]))
	abs_acc = (u[0]+u[1]+u[2]+u[3])/mass
	vel = vertcat(x[3],x[4],x[5])
	ang_x = x[6]
	ang_y = x[7]
	ang_z = x[8]
	angrate_x = x[9]
	angrate_y = x[10]
	angrate_z = x[11]
	cp = cos(-ang_y)
	sp = sin(-ang_y)
	sr = sin(-ang_x)
	cr = cos(-ang_x)
	sy = sin(ang_z)
	cy = cos(ang_z)
	R_B2W_11 = cp * cy
	R_B2W_21 = ((sr * sp * cy) - (cr * sy))
	R_B2W_31 = ((cr * sp * cy) + (sr * sy))
	R_B2W_12 = cp * sy
	R_B2W_22 = ((sr * sp * sy) + (cr * cy))
	R_B2W_32 = ((cr * sp * sy) - (sr * cy))
	R_B2W_13 = -sp
	R_B2W_23 = sr * cp
	R_B2W_33 = cr * cp
	new_acc_x = abs_acc * R_B2W_13
	new_acc_y = abs_acc * R_B2W_23
	new_acc_z = abs_acc * R_B2W_33 - gravity
	new_acc = vertcat(new_acc_x, new_acc_y, new_acc_z)
	angrate_W_x = R_B2W_11*angrate_x + R_B2W_12*angrate_y + R_B2W_13*angrate_z
	angrate_W_y = R_B2W_21*angrate_x + R_B2W_22*angrate_y + R_B2W_23*angrate_z
	angrate_W_z = R_B2W_31*angrate_x + R_B2W_32*angrate_y + R_B2W_33*angrate_z
	angrate_W=vertcat(angrate_W_x,angrate_W_y,angrate_W_z)
	xdot = vertcat(vel, new_acc, angrate_W, dangvel_Body_dt)*t_trans
	return xdot
def rk4(ode, h, t, x, u):
	k1 = ode(t,x,u);
	k2 = ode(t,x+h/2*k1,u);
	k3 = ode(t,x+h/2*k2,u);
	k4 = ode(t,x+h*k3,  u);
	xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
	return xf
def handle_traj_gen(req):
	vel_constr = 2
	angle_constr = 3.14159265/4
	yaw_init_state = 0
	yaw_locked = True
	if yaw_locked:
		yaw_constr_plus = yaw_init_state
		yaw_constr_minus = yaw_init_state
	else:
		yaw_constr_plus = 3.14159265*4
		yaw_constr_minus = -3.14159265*4
	input_constr = 6
	x0pos=[0,0,1]
	x0vel=[0,0,0]
	x0ang=[0,0,yaw_init_state]
	x0angrate=[0,0,0]
	x1pos = [2.,0.,2.5]
	x1vel_y = 0
	x1vel_z = 0
	x2pos = [2.,4.,2.5]
#??? x y z mistake
	x2vel_y = 0
	x2vel_z = 0
	xFpos = [0,0,1.2]
	xFvel = [0,0,0]
	xFang = [0,0,yaw_init_state]
	xFangrate=[0,0,0]
	N = 25
# x = vertcat(pos,vel,ang,angrate)
	x=MX.sym('x',12)
	u=MX.sym('u',4)
	t_transf = MX.sym('t_transf')
	t = MX.sym('t')
	
	xdot = diff_eq_time_trans(0,x,u,t_transf)
	L = (u[0]**2 + u[1]**2 + u[2]**2 + u[3]**2)*t_transf

	f = Function('f', [t, x, u, t_transf],[xdot])
	fq = Function('fq', [t, x, u, t_transf],[L])

	X0 = MX.sym('X0', 12)
	U=MX.sym('U',4)
	TF = MX.sym('TF')
	M = 1
	DT = 1. / N / M
	X=X0
	Q=0
	f_handle = lambda t,x,u: f(t,x,u,TF)
	fq_handle = lambda t,x,u: fq(t,x,u,TF)
	for j in range(M):
		X=rk4(f_handle, DT, 0., X, U)
		Q=rk4(fq_handle, DT, 0., 0, U)
	F = Function('F', [X0, U, TF], [X, Q],['x0','p','tf_interval'],['xf','qf'])
	# Start with an empty NLP
	w=[]
	w0 = []
	lbw = []
	ubw = []
	J = 0
	g=[]
	lbg = []
	ubg = []
	# "Lift" initial conditions
#
	t_interval_1 = MX.sym('t_interval_1')
	t_interval_2 = MX.sym('t_interval_2')
	t_interval_3 = MX.sym('t_interval_3')
#

	Xk = MX.sym('X0', 12)
	w += [Xk]
	# this initial condition has the specified bound, 
	# with only one possible solution
	# x0(0)=1; x1(0)=0;
	lbw += x0pos + x0vel + x0ang + x0angrate
	ubw += x0pos + x0vel + x0ang + x0angrate
	#solver's guess
	w0 += x0pos + x0vel + x0ang + x0angrate
	gravity, mass, kF, kM, Len = project_parameters()
	for k in range(N):
		Uk = MX.sym('U_' + str(k), 4)
		w   += [Uk]
		if k!=0:
			lbw += [0,0,0,0]
			ubw += [input_constr,input_constr,input_constr,input_constr]
		else:
			lbw += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]
			ubw += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]
		w0  += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]

		Fk = F(x0=Xk, p=Uk, tf_interval=t_interval_1)
		Xk_end = Fk['xf']
		J=J+Fk['qf']

		Xk = MX.sym('X_' + str(k+1), 12)
		w   += [Xk]
		if k!=N-1:
			lbw += [-inf,-inf,-inf, -vel_constr,-vel_constr,-vel_constr, -angle_constr,-angle_constr,yaw_constr_minus, -inf,-inf,-inf]
			ubw += [inf,inf,inf, vel_constr,vel_constr,vel_constr, angle_constr,angle_constr,yaw_constr_plus, inf,inf,inf]
		else:
			lbw += x1pos + [-vel_constr,x1vel_y,x1vel_z, -angle_constr,-angle_constr,yaw_constr_minus, -inf,-inf,-inf]
			ubw += x1pos +  [vel_constr,x1vel_y,x1vel_z, angle_constr,angle_constr,yaw_constr_plus, inf,inf,inf]
		w0  += x0pos+[0,0,0, 0,0,3, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]

	for k in range(N,2*N):
		Uk = MX.sym('U_' + str(k), 4)
		w   += [Uk]
		lbw += [0,0,0,0]
		ubw += [input_constr,input_constr,input_constr,input_constr]
		w0  += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]

		Fk = F(x0=Xk, p=Uk, tf_interval=t_interval_2)
		Xk_end = Fk['xf']
		J=J+Fk['qf']

		Xk = MX.sym('X_' + str(k+1), 12)
		w   += [Xk]
		if k!=2*N-1:
			lbw += [-inf,-inf,-inf, -vel_constr,-vel_constr,-vel_constr, -angle_constr,-angle_constr,yaw_constr_minus, -inf,-inf,-inf]
			ubw += [inf,inf,inf, vel_constr,vel_constr,vel_constr, angle_constr,angle_constr,yaw_constr_plus, inf,inf,inf]
		else:
#??? x y z mistake
			lbw += x2pos + [-vel_constr,x2vel_y,x2vel_z, -angle_constr,-angle_constr,yaw_constr_minus, -inf,-inf,-inf]
			ubw += x2pos + [vel_constr,x2vel_y,x2vel_z, angle_constr,angle_constr,yaw_constr_plus, inf,inf,inf]
		w0  += x1pos+[0,0,0, 0,0,3, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
	for k in range(2*N,3*N):
		Uk = MX.sym('U_' + str(k), 4)
		w   += [Uk]
		if k!=3*N-1:
			lbw += [0,0,0,0]
			ubw += [input_constr,input_constr,input_constr,input_constr]
		else:
			lbw += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]
			ubw += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]
		w0  += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]

		Fk = F(x0=Xk, p=Uk, tf_interval=t_interval_3)
		Xk_end = Fk['xf']
		J=J+Fk['qf']

		Xk = MX.sym('X_' + str(k+1), 12)
		w   += [Xk]
		if k!=3*N-1:
			lbw += [-inf,-inf,-inf, -vel_constr,-vel_constr,-vel_constr, -angle_constr,-angle_constr,yaw_constr_minus, -inf,-inf,-inf]
			ubw += [inf,inf,inf, vel_constr,vel_constr,vel_constr, angle_constr,angle_constr,yaw_constr_plus, inf,inf,inf]
		else:
			lbw += xFpos + xFvel + xFang + xFangrate
			ubw += xFpos + xFvel + xFang + xFangrate
		w0  += x2pos + [0,0,0, 0,0,3, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
	w += [t_interval_1,t_interval_2,t_interval_3]
	lbw += [0.1,0.1,0.1]
	ubw += [inf,inf,inf]
	w0 += [2,3,2]
	# g += [t_interval_1,t_interval_2,t_interval_3]
	# lbg += [0,0,0]
	# ubg += [inf,inf,inf]
	prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
	solver = nlpsol('solver', 'ipopt', prob)
	sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
	w_opt_all = sol['x'].full().flatten()

	t1_opt = w_opt_all[-3]
	t2_opt = w_opt_all[-2]
	t3_opt = w_opt_all[-1]
	print(t1_opt)
	print(t2_opt)
	print(t3_opt)
	w_opt = np.delete(w_opt_all, [-3,-2,-1])
	pos_x_opt = w_opt[0::16]
	pos_y_opt = w_opt[1::16]
	pos_z_opt = w_opt[2::16]
	vel_x_opt = w_opt[3::16]
	vel_y_opt = w_opt[4::16]
	vel_z_opt = w_opt[5::16]
	ang_x_opt = w_opt[6::16]
	ang_y_opt = w_opt[7::16]
	ang_z_opt = w_opt[8::16]
	angr_x_opt = w_opt[9::16]
	angr_y_opt = w_opt[10::16]
	angr_z_opt = w_opt[11::16]
	u1_opt=[]
	u2_opt=[]
	u3_opt=[]
	u4_opt=[]
	for k in range(3*N):
		u1_opt += [w_opt[12+k*16]]
		u2_opt += [w_opt[13+k*16]]
		u3_opt += [w_opt[14+k*16]]
		u4_opt += [w_opt[15+k*16]]
	tgrid1 = [t1_opt/N*k for k in range(N)]
	tgrid2 = [t1_opt+t2_opt/N*k for k in range(N)]
	tgrid3 = [t1_opt+t2_opt+t3_opt/N*k for k in range(N+1)]
	tgrid3u = [t1_opt+t2_opt+t3_opt/N*k for k in range(N)]
	tgrid = tgrid1+tgrid2+tgrid3
	tgridu = tgrid1+tgrid2+tgrid3u

	msg=MultiDOFJointTrajectory()
	msg.header.frame_id =''
	msg.header.stamp = rospy.Time.now()
	msg.joint_names = ["base_link"]
	for i in range(len(tgrid)):
		transforms = Transform()
		transforms.translation.x = pos_x_opt[i]
		transforms.translation.y = pos_y_opt[i]
		transforms.translation.z = pos_z_opt[i]

		quaternion = tf.transformations.quaternion_from_euler(ang_x_opt[i], ang_y_opt[i], ang_z_opt[i])
		transforms.rotation = Quaternion(quaternion[0],quaternion[1],quaternion[2],quaternion[3])

		velocities = Twist()
		velocities.linear.x = vel_x_opt[i]
		velocities.linear.y = vel_y_opt[i]
		velocities.linear.z = vel_z_opt[i]
		# velocities.angular.x = 
		accelerations=Twist()
		# time_from_start += rospy.Duration.from_sec(waypoints[i].waiting_time)
		point = MultiDOFJointTrajectoryPoint([transforms],[velocities],[accelerations],rospy.Duration.from_sec(tgrid[i]))
		
		msg.points.append(point)
	return msg

def traj_gen_server():
	rospy.init_node('traj_gen_server')
	s = rospy.Service('traj_gen', traj_gen_service, handle_traj_gen)
	print "Ready to generate trajectory"
	rospy.spin()
if __name__ == "__main__":
	traj_gen_server()
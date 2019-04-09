from casadi import *
import casadi.tools as ctools
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from geometry_msgs.msg import Transform, Quaternion
import tf
def project_parameters():
	kF=10
	kM=2
	mass=1.0
	Len=0.225
	gravity=9.81
	return gravity, mass, kF, kM, Len

if __name__=='__main__':
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
	input_constr = 4
	x0pos=[0,0,1]
	x0vel=[0,0,0]
	x0ang=[0,0,yaw_init_state]
	x0angrate=[0,0,0]
	x1pos = [2.,0.,2.5]
	x1vel_y = 0
	x1vel_z = 0
	x2pos = [6.,4.,2.5]
	x2vel_x = 0
	x2vel_z = 0
	xFpos = [7.,6.,1.2]
	xFvel = [0,0,0]
	xFang = [0,0,yaw_init_state]
	xFangrate=[0,0,0]
	N = 25
	pos_x = MX.sym('pos_x')
	pos_y = MX.sym('pos_y')
	pos_z = MX.sym('pos_z')
	pos=vertcat(pos_x,pos_y,pos_z)
	vel_x = MX.sym('vel_x')
	vel_y = MX.sym('vel_y')
	vel_z = MX.sym('vel_z')
	vel=vertcat(vel_x,vel_y,vel_z)
	ang_x = MX.sym('ang_x')
	ang_y = MX.sym('ang_y')
	ang_z = MX.sym('ang_z')
	ang=vertcat(ang_x,ang_y,ang_z)
	angrate_x = MX.sym('angrate_x')
	angrate_y = MX.sym('angrate_y')
	angrate_z = MX.sym('angrate_z')
	angrate=vertcat(angrate_x,angrate_y,angrate_z)
	x = vertcat(pos,vel,ang,angrate)

	u1 = MX.sym('u1')
	u2 = MX.sym('u2')
	u3 = MX.sym('u3')
	u4 = MX.sym('u4')
	u = vertcat(u1, u2, u3, u4)
#
	t_transf = MX.sym('t_transf')
#
	gravity, mass, kF, kM, Len = project_parameters()
	dangvel_Body_dt = vertcat(kF*Len*(u2-u4), kF*Len*(-u1+u3), kM*(u1-u2+u3-u4))
	abs_acc = (u1+u2+u3+u4)/mass
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
	# R_B2W_13 = ((cr * sp * cy) + (sr * sy))
	# R_B2W_23 = ((cr * sp * sy) - (sr * cy))
	# R_B2W_33 = cr * cp
	new_acc_x = abs_acc * R_B2W_13
	new_acc_y = abs_acc * R_B2W_23
	new_acc_z = abs_acc * R_B2W_33 - gravity
	new_acc = vertcat(new_acc_x, new_acc_y, new_acc_z)
	angrate_W_x = R_B2W_11*angrate_x + R_B2W_12*angrate_y + R_B2W_13*angrate_z
	angrate_W_y = R_B2W_21*angrate_x + R_B2W_22*angrate_y + R_B2W_23*angrate_z
	angrate_W_z = R_B2W_31*angrate_x + R_B2W_32*angrate_y + R_B2W_33*angrate_z
	angrate_W=vertcat(angrate_W_x,angrate_W_y,angrate_W_z)
	xdot = vertcat(vel, new_acc, angrate_W, dangvel_Body_dt)*t_transf
#
	L = (u1**2 + u2**2 + u3**2 + u4**2)*t_transf
#
	M = 1
	
	f = Function('f', [x, u, t_transf],[xdot, L])
	X0 = MX.sym('X0', 12)
	U=MX.sym('U',4)
	TF = MX.sym('TF')
#
	DT = 1. / N / M
#
	X=X0
	Q=0
	for j in range(M):
		k1, k1_q = f(X, U, TF)
		k2, k2_q = f(X + DT/2 * k1, U, TF)
		k3, k3_q = f(X + DT/2 * k2, U, TF)
		k4, k4_q = f(X + DT * k3, U, TF)
		X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
		Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
	F = Function('F', [X0, U, TF], [X, Q],['x0','p','tf_interval'],['xf','qf'])

	# Evaluate at a test point
	# Fk = F(x0=[0.2,0.3],p=0.4)
	# print(Fk['xf'])
	# print(Fk['qf'])

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

	for k in range(N):
		Uk = MX.sym('U_' + str(k), 4)
		w   += [Uk]
		lbw += [0,0,0,0]
		ubw += [input_constr,input_constr,input_constr,input_constr]
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
			lbw += x2pos + [x2vel_x,-vel_constr,x2vel_z, -angle_constr,-angle_constr,yaw_constr_minus, -inf,-inf,-inf]
			ubw += x2pos + [x2vel_x,vel_constr,x2vel_z, angle_constr,angle_constr,yaw_constr_plus, inf,inf,inf]
		w0  += x1pos+[0,0,0, 0,0,3, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0, 0,0,0]
	for k in range(2*N,3*N):
		Uk = MX.sym('U_' + str(k), 4)
		w   += [Uk]
		lbw += [0,0,0,0]
		ubw += [input_constr,input_constr,input_constr,input_constr]
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
	quat_opt_x = []
	quat_opt_y = []
	quat_opt_z = []
	quat_opt_w = []
	for i in range(len(ang_x_opt)):
		quaternion = tf.transformations.quaternion_from_euler(ang_x_opt[i], ang_y_opt[i], ang_z_opt[i])
		quat_opt_x += [quaternion[0]]
		quat_opt_y += [quaternion[1]]
		quat_opt_z += [quaternion[2]]
		quat_opt_w += [quaternion[3]]

	angr_x_opt = w_opt[9::16]
	angr_y_opt = w_opt[10::16]
	angr_z_opt = w_opt[11::16]
	u1_opt=[]
	u2_opt=[]
	u3_opt=[]
	u4_opt=[]
	for k in range(3*N):
		u1_opt += [w_opt[12+k*16]]
		u2_opt += [-w_opt[13+k*16]]
		u3_opt += [w_opt[14+k*16]]
		u4_opt += [-w_opt[15+k*16]]
	tgrid1 = [t1_opt/N*k for k in range(N)]
	tgrid2 = [t1_opt+t2_opt/N*k for k in range(N)]
	tgrid3 = [t1_opt+t2_opt+t3_opt/N*k for k in range(N+1)]
	tgrid3u = [t1_opt+t2_opt+t3_opt/N*k for k in range(N)]
	tgrid = tgrid1+tgrid2+tgrid3
	tgridu = tgrid1+tgrid2+tgrid3u

	fig = plt.figure(1)
	plt.clf()
	pos_subplot = fig.add_subplot(511)
	pos_subplot.plot(tgrid, pos_x_opt, '--')
	pos_subplot.plot(tgrid, pos_y_opt, '--')
	pos_subplot.plot(tgrid, pos_z_opt, '--')
	pos_subplot.set_xlabel('t')
	pos_subplot.legend(['x','y','z'])

	vel_subplot = fig.add_subplot(512)
	vel_subplot.plot(tgrid, vel_x_opt, '-')
	vel_subplot.plot(tgrid, vel_y_opt, '-')
	vel_subplot.plot(tgrid, vel_z_opt, '-')
	vel_subplot.set_xlabel('t')
	vel_subplot.legend(['vx','vy','vz'])
	ang_subplot = fig.add_subplot(513)
	ang_subplot.plot(tgrid, quat_opt_x, '-')
	ang_subplot.plot(tgrid, quat_opt_y, '-')
	ang_subplot.plot(tgrid, quat_opt_z, '-')
	ang_subplot.set_xlabel('t')
	ang_subplot.legend(['roll','pitch','yaw'])
	angr_subplot = fig.add_subplot(514)
	angr_subplot.plot(tgrid, angr_x_opt, '-')
	angr_subplot.plot(tgrid, angr_y_opt, '-')
	angr_subplot.plot(tgrid, angr_z_opt, '-')
	angr_subplot.set_xlabel('t')
	angr_subplot.legend(['roll','pitch','yaw'])
	rotor_subplot = fig.add_subplot(515)
	rotor_subplot.plot(tgridu, u1_opt,'-')
	rotor_subplot.plot(tgridu, u2_opt,'-')
	rotor_subplot.plot(tgridu, u3_opt,'-')
	rotor_subplot.plot(tgridu, u4_opt,'-')
	rotor_subplot.set_xlabel('t')
	rotor_subplot.legend(['1','2','3','4'])
	# rotor_subplot.xlabel('t')
	# rotor_subplot.legend(['1','2','3','4'])
	plt.grid()
	plt.show()
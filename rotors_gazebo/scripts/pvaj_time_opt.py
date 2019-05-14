#!/usr/bin/env python
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
def diff_eq_time_trans(tau,x,u,t_trans):
	vel = vertcat(x[3],x[4],x[5])
	acc = vertcat(x[6],x[7],x[8])
	dacc_dt = vertcat(u[0],u[1],u[2])

	xdot = vertcat(vel,acc,dacc_dt) * t_trans
	return xdot
def rk4(ode, h, t, x, u):
	k1 = ode(t,x,u);
	k2 = ode(t,x+h/2*k1,u);
	k3 = ode(t,x+h/2*k2,u);
	k4 = ode(t,x+h*k3,  u);
	xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
	return xf
def plot_rings(radius, yaw, x_tr, y_tr, z_tr, figure_num):
	pass
if __name__=='__main__':
	vel_constr = 3
	acc_constr = 5
	input_constr = 10

	x0pos=[0,0,0]
	x0vel=[0,0,0]
	x0acc=[0,0,0]
	x1pos = [2.,0.,1.5]
	x1vel_y = 0
	x1vel_z = 0
	x2pos = [6.,4.,1.5]
	x2vel_x = 0
	x2vel_z = 0
	xFpos = [7,6,0.3]
	xFvel = [0,0,0]
	xFacc = [0,0,0]

	N = 25
	x=MX.sym('x',9)
	u=MX.sym('u',3)
	t_transf = MX.sym('t_transf')
	t = MX.sym('t')
	xdot = diff_eq_time_trans(0,x,u,t_transf)
	Lagr_term = (u[0]**2 + u[1]**2 + u[2]**2) * t_transf
	f = Function('f', [t, x, u, t_transf],[xdot])
	fq = Function('fq', [t, x, u, t_transf],[Lagr_term])

	X0 = MX.sym('X0', 9)
	U=MX.sym('U',3)
	TF = MX.sym('TF')
	# substeps
	M = 1
	# step length for rk4
	DT = 1. / N / M
	X=X0
	Q=0
	f_handle = lambda t,x,u: f(t,x,u,TF)
	fq_handle = lambda t,x,u: fq(t,x,u,TF)
	for j in range(M):
		#rk4(ode, h, t, x, u)
		X=rk4(f_handle, DT, 0., X, U)
		Q=rk4(fq_handle, DT, 0., 0, U)
	# function that combines the computation of diff eq and Lagr term
	F = Function('F', [X0, U, TF], [X, Q],['x0','p','tf_interval'],['xf','qf'])
	# Start with an empty NLP
	# w contains all u and x
	w=[]
	#guess for solver
	w0 = []
	lbw = []
	ubw = []
	J = 0
	# g gives equational constrains such as Xk_end-Xk == 0
	g=[]
	lbg = []
	ubg = []

	t_interval_1 = MX.sym('t_interval_1')
	t_interval_2 = MX.sym('t_interval_2')
	t_interval_3 = MX.sym('t_interval_3')

	Xk = MX.sym('X0', 9)
	w += [Xk]
	# this initial condition has the specified bound, 
	# with only one possible solution
	# x0(0)=1; x1(0)=0;
	lbw += x0pos + x0vel + x0acc
	ubw += x0pos + x0vel + x0acc
	#solver's guess
	w0 += x0pos + x0vel + x0acc

	for k in range(N):
		Uk = MX.sym('U_' + str(k), 3)
		w   += [Uk]
		if k!=0:
			lbw += [-input_constr,-input_constr,-input_constr]
			ubw += [input_constr,input_constr,input_constr]
		else:
			lbw += [0,0,0]
			ubw += [0,0,0]
		w0  += [0,0,0]
		# function that combines the computation of diff eq and Lagr term
		Fk = F(x0=Xk, p=Uk, tf_interval=t_interval_1)
		Xk_end = Fk['xf']
		J=J+Fk['qf']

		Xk = MX.sym('X_' + str(k+1), 9)
		w   += [Xk]
		if k!=N-1:
			lbw += [-inf,-inf,-inf, -vel_constr,-vel_constr,-vel_constr, -acc_constr,-acc_constr,-acc_constr]
			ubw += [inf,inf,inf, vel_constr,vel_constr,vel_constr, acc_constr,acc_constr,acc_constr]
		else:
			lbw += x1pos + [-vel_constr,x1vel_y,x1vel_z, -acc_constr,-acc_constr,-acc_constr]
			ubw += x1pos +  [vel_constr,x1vel_y,x1vel_z, acc_constr,acc_constr,acc_constr]
		w0  += x0pos+[0,0,0, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0]
	for k in range(N,2*N):
		Uk = MX.sym('U_' + str(k), 3)
		w   += [Uk]
	#	if k!=0:
		lbw += [-input_constr,-input_constr,-input_constr]
		ubw += [input_constr,input_constr,input_constr]
	#	else:
	#		lbw += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]
	#		ubw += [0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass,0.25*gravity*mass]
		w0  += [0,0,0]
		# function that combines the computation of diff eq and Lagr term
		Fk = F(x0=Xk, p=Uk, tf_interval=t_interval_2)
		Xk_end = Fk['xf']
		J=J+Fk['qf']

		Xk = MX.sym('X_' + str(k+1), 9)
		w   += [Xk]
		if k!=2*N-1:
			lbw += [-inf,-inf,-inf, -vel_constr,-vel_constr,-vel_constr, -acc_constr,-acc_constr,-acc_constr]
			ubw += [inf,inf,inf, vel_constr,vel_constr,vel_constr, acc_constr,acc_constr,acc_constr]
		else:
			lbw += x2pos + [x2vel_x,-vel_constr,x2vel_z, -acc_constr,-acc_constr,-acc_constr]
			ubw += x2pos +  [x2vel_x,vel_constr,x2vel_z, acc_constr,acc_constr,acc_constr]
		w0  += x1pos+[0,0,0, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0]
	for k in range(2*N,3*N):
		Uk = MX.sym('U_' + str(k), 3)
		w   += [Uk]
		if k!=3*N-1:
			lbw += [-input_constr,-input_constr,-input_constr]
			ubw += [input_constr,input_constr,input_constr]
		else:
			lbw += [0,0,0]
			ubw += [0,0,0]
		w0  += [0,0,0]
		# function that combines the computation of diff eq and Lagr term
		Fk = F(x0=Xk, p=Uk, tf_interval=t_interval_3)
		Xk_end = Fk['xf']
		J=J+Fk['qf']

		Xk = MX.sym('X_' + str(k+1), 9)
		w   += [Xk]
		if k!=3*N-1:
			lbw += [-inf,-inf,-inf, -vel_constr,-vel_constr,-vel_constr, -acc_constr,-acc_constr,-acc_constr]
			ubw += [inf,inf,inf, vel_constr,vel_constr,vel_constr, acc_constr,acc_constr,acc_constr]
		else:
			lbw += xFpos + xFvel + xFacc
			ubw += xFpos + xFvel + xFacc
		w0  += x2pos+[0,0,0, 0,0,0]
		g   += [Xk_end-Xk]
		lbg += [0,0,0, 0,0,0, 0,0,0]
		ubg += [0,0,0, 0,0,0, 0,0,0]
	w += [t_interval_1,t_interval_2,t_interval_3]
	lbw += [0.1,0.1,0.1]
	ubw += [inf,inf,inf]
	w0 += [2,3,2]
	J_time = t_interval_1+t_interval_2+t_interval_3 + J/1000
	prob = {'f': J_time, 'x': vertcat(*w), 'g': vertcat(*g)}
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
	pos_x_opt = w_opt[0::12]
	pos_y_opt = w_opt[1::12]
	pos_z_opt = w_opt[2::12]
	vel_x_opt = w_opt[3::12]
	vel_y_opt = w_opt[4::12]
	vel_z_opt = w_opt[5::12]
	acc_x_opt = w_opt[6::12]
	acc_y_opt = w_opt[7::12]
	acc_z_opt = w_opt[8::12]
	jerk_x_opt=[]
	jerk_y_opt=[]
	jerk_z_opt=[]
	for k in range(3*N):
		jerk_x_opt += [w_opt[9+k*12]]
		jerk_y_opt += [-w_opt[10+k*12]]
		jerk_z_opt += [w_opt[11+k*12]]
#	jerk_x_opt = w_opt[9::12]
#	jerk_y_opt = w_opt[10::12]
#	jerk_z_opt = w_opt[11::12]

	tgrid1 = [t1_opt/N*k for k in range(N)]
	tgrid2 = [t1_opt+t2_opt/N*k for k in range(N)]
	tgrid3 = [t1_opt+t2_opt+t3_opt/N*k for k in range(N+1)]
	tgrid3u = [t1_opt+t2_opt+t3_opt/N*k for k in range(N)]
	tgrid = tgrid1+tgrid2+tgrid3
	tgridu = tgrid1+tgrid2+tgrid3u

	fig = plt.figure(1)
	plt.clf()
	pos_subplot = fig.add_subplot(411)
	pos_subplot.plot(tgrid, pos_x_opt, '--')
	pos_subplot.plot(tgrid, pos_y_opt, '--')
	pos_subplot.plot(tgrid, pos_z_opt, '--')
	pos_subplot.set_xlabel('t')
	pos_subplot.legend(['x','y','z'])

	vel_subplot = fig.add_subplot(412)
	vel_subplot.plot(tgrid, vel_x_opt, '-')
	vel_subplot.plot(tgrid, vel_y_opt, '-')
	vel_subplot.plot(tgrid, vel_z_opt, '-')
	vel_subplot.set_xlabel('t')
	vel_subplot.legend(['vx','vy','vz'])
	ang_subplot = fig.add_subplot(413)
	ang_subplot.plot(tgrid, acc_x_opt, '-')
	ang_subplot.plot(tgrid, acc_y_opt, '-')
	ang_subplot.plot(tgrid, acc_z_opt, '-')
	ang_subplot.set_xlabel('t')
	ang_subplot.legend(['ax','ay','az'])
	angr_subplot = fig.add_subplot(414)
	angr_subplot.plot(tgridu, jerk_x_opt, '-')
	angr_subplot.plot(tgridu, jerk_y_opt, '-')
	angr_subplot.plot(tgridu, jerk_z_opt, '-')
	angr_subplot.set_xlabel('t')
	angr_subplot.legend(['jx','jy','jz'])
	plt.grid()
	plt.show()
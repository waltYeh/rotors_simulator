function uav_time_opt()
addpath('../../../casadi')
import casadi.*;

clear;
clc;
x0=zeros(9,1);
x0(7:9)=[0;0;3];
x1pos=[2;0;1.5];
x1vel=[nan;0;0];
x2pos=[6;4;1.5];
x2vel=[0;nan;0];
xF=zeros(9,1);
xF(1:3)=[7;6;1];
xF(7:9)=[0;0;3];

% Q=diag([1,1,1,1]);
R=diag([1e-6,1e-6,1e-6]);
% Q=diag([100,100,10,10]);
% R=diag([1000,10]);

% [~, ~, ~, ~, state_constr, input_constr] = project_parameters;
vel_constr = 5;
acc_constr = 5;
input_constr = 10;
N=100;
opti=casadi.Opti();
X=opti.variable(9,3*N+1);
U=opti.variable(3,3*N);
t_interval=opti.variable(3,1);
tau=0;
J1=0;
for i=1:N
    tau=tau+1/N;
%     x_next = rk4(@(t,x,u)diff_eq(t,x,u),tf/N,tau,X(:,i),U(:,i));
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,t_interval(1)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J1 = J1 + U(:,i)' * R * U(:,i)*1/N;
end
opti.subject_to(X(:,1) == x0);
opti.subject_to(X(1:3,N+1) == x1pos);
opti.subject_to(X(5:6,N+1) == x1vel(2:3));
J2=0;
for i=N+1:2*N
    tau=tau+1/N;
%     x_next = rk4(@(t,x,u)diff_eq(t,x,u),tf/N,tau,X(:,i),U(:,i));
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,t_interval(2)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J2 = J2 + U(:,i)' * R * U(:,i)*1/N;
end
% opti.subject_to(X(:,1) == x0);
opti.subject_to(X(1:3,2*N+1) == x2pos);
opti.subject_to(X(4,2*N+1) == x2vel(1));
opti.subject_to(X(6,2*N+1) == x2vel(3));

J3=0;
for i=2*N+1:3*N
    tau=tau+1/N;
%     x_next = rk4(@(t,x,u)diff_eq(t,x,u),tf/N,tau,X(:,i),U(:,i));
    x_next = rk4(@(t,x,u)diff_eq_time_trans(t,x,u,t_interval(3)),1/N,tau,X(:,i),U(:,i));
    opti.subject_to(X(:,i+1)==x_next);
    J3 = J3 + U(:,i)' * R * U(:,i)*1/N;
end

opti.subject_to(X(:,end) == xF);
opti.subject_to(-vel_constr <= X(4:6,:) <= vel_constr);
opti.subject_to(-acc_constr <= X(7:9,:) <= acc_constr);
opti.subject_to(-input_constr <= U(1:3,:) <= input_constr);
opti.subject_to(t_interval(1)> 0);
opti.subject_to(t_interval(2)> 0);
opti.subject_to(t_interval(3)> 0);
opti.minimize(t_interval(1)+t_interval(2)+t_interval(3)+J1+J2+J3);
opti.solver('ipopt');
opti.set_initial(t_interval,[3;3;3]);
sol = opti.solve();
pos=sol.value(X(1:3,:))';
vel=sol.value(X(4:6,:))';
acc=sol.value(X(7:9,:))';
jerk=sol.value(U(1:3,:))';
tx1=linspace(0,sol.value(t_interval(1)),N+1);
tx2=linspace(0,sol.value(t_interval(2)),N+1)+tx1(end);
tx3=linspace(0,sol.value(t_interval(3)),N+1)+tx2(end);
tx=[tx1,tx2(2:end),tx3(2:end)];
tu=tx(1:end-1);
figure(1)
subplot(4,1,1)
plot(tx,pos)
title('pos')
xlabel('t (s)')
ylabel('pos (m)')
legend('x','y','z')
grid
subplot(4,1,2)
plot(tx,vel)
title('vel')
xlabel('t (s)')
ylabel('velocity (m/s)')
legend('x','y','z')
grid
subplot(4,1,3)
plot(tx,acc)
title('acc')
xlabel('t (s)')
ylabel('acc (m/s^2)')
legend('x','y','z')
grid
subplot(4,1,4)
plot(tu,jerk)
title('jerk')
xlabel('t (s)')
ylabel('jerk (m/s^3)')
legend('x','y','z')
grid
figure(2)
plot3(pos(:,1),pos(:,2),pos(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
grid
end
function xf = rk4(ode,h,t,x,u)
  k1 = ode(t,x,u);
  k2 = ode(t,x+h/2*k1,u);
  k3 = ode(t,x+h/2*k2,u);
  k4 = ode(t,x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end
function dx=diff_eq_time_trans(tau,x,u,tf)
pos=x(1:3);
vel=x(4:6);
acc=x(7:9);

dacc_dt=u;
dvel_dt=acc;
dpos_dt=vel;

dx=tf*[dpos_dt;dvel_dt;dacc_dt];
end
function dx=diff_eq(t,x,u)
pos=x(1:3);
vel=x(4:6);
acc=x(7:9);

dacc_dt=u;
dvel_dt=acc;
dpos_dt=vel;

dx=[dpos_dt;dvel_dt;dacc_dt];
end
function [state_constr, input_constr] = project_parameters
%% Definition of system parameters


%% Constraints
state_constr=3/2*pi;
input_constr=1000;

%% Scaling factors

end


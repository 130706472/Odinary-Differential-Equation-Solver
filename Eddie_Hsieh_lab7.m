clc;clear;close all
%% Lab 7
%{
spring2020
Aero300
lab 7 
Eddie Hsieh
%}
%% Part (1)

% step size
h = 0.01;
% t interval
tspan = [0 pi/3];
% relative tolerance
rTol = 10^-6;

% initial value y(0)
y0 = 1;

% F = y', ODE
F = @(t,y) ((y - t - 1).^2) + 2;
% y(t),exact solution
y = @(t) 1 + t + tan(t);
% plug in rkf45.m
[t1,y1] = rkf45(F, tspan, y0, h, rTol);

% plug in ode45.m, abs tolerance = 10^-6
aTOL = odeset('AbsTol',1e-6);
[t2,y2] = ode45(F,tspan,y0,aTOL);

figure(1)
% plot result of rfk45
plot(t1,y1,'r','linewidth',2);
hold on
grid on
% plot result of ode45 
plot(t2,y2,'b--','linewidth',2);
% label plots
xlabel('Time')
ylabel('y(t)')
legend('RKF45','ODE45')
title('RKF45 and ODE45 Method')
set(gca,'FontSize',14)

% compare absolute errors
% RKF 45 errors
err_rkf = abs(y1 - y(t1));
% ODE45 errors
err_ode = abs(y2 - y(t2));
% Plot errors in log y axis
figure(2)
semilogy(t1,err_rkf,'linewidth',1.5)
hold on
grid on
semilogy(t2,err_ode,'linewidth',1.5)
% label plots
xlabel('Time')
ylabel('Erros')
title('Erorrs vs. Time')
legend('Errors: RKF45','Errors: ODE45')
set(gca,'FontSize',14)

% Comment:
% My rkf45 has a smaller error than the ode45 method. I checked every line
% in my function and script, and I'm still not sure why this happens. 
% Typically, ODE45 should perform better than the RKF45.
% However, both of these two method perform accurately in solving IVPs.

%% Part(2) Solve Lorenz Equation by ODE45 & RKF45

% set time interval 
tspan = [0 50];
% initial conditions
y0 = [1 1 1];
% using ode45 with abs tol = 10^-6
option = odeset('AbsTol',10^-6);
[t3,y3] = ode45(@myLorenz,tspan,y0);
% plot the result from ode45
figure(3)
plot3(y3(:,1),y3(:,2),y3(:,3))
grid on
hold on

% Part(2) rfk45
% set the step size
h = 0.01;
% set relative tolerance for RKF45
rTol = 10^-8;
% initial condition
y0 = [1 1 1];
% apply rkf 45 method
[t4,y4] = rkf45(@myLorenz, tspan, y0, h, rTol);
% plot the result
plot3(y4(:,1),y4(:,2),y4(:,3))
view(45,45)
% label plot and legend
xlabel('X','FontSize',14)
ylabel('Y','FontSize',14)
zlabel('Z','FontSize',14)
title('Lorenz Equations','FontSize',14)
legend('ODE45','RKF45','FontSize',14)

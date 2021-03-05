%% Lorenz differential Equation
function F = myLorenz(t,x)
%{
Eddie Hsieh
Lorenz equation function

Input: 
    x: vector of [x y z] values
    t: time( not a variable in this function because of non-linear system)
Output:
    F: vector of [dx/dt dy/dt dz/dt]
%}

% set parameter
% rho sigma beta
rho = 28;
sigma = 10;
beta = 8/3;

% differential equations

dxdt = sigma*(x(2) - x(1));

dydt = x(1)*(rho - x(3)) - x(2); 

dzdt = x(1)*x(2) - beta*x(3);

% F = [dxdt dydt dzdt];
F = [dxdt;dydt;dzdt];
end
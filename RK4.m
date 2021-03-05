function [t, w] = RK4(f, t0, y0, T, n)

% Implements the 4th-order Runge-Kutta method for 1st-order IVPs
%   f: right-hand side of ODE
%   t0: Initial time
%   y0: Function value at the initial time
%   T: final time at which estimate is desired
%   n: number of iterations needed

h = (T-t0)/n;
h2 = h/2;
t = linspace(t0, T, n+1)';

w = zeros(n+1,1);
w(1) = y0;

for i=1:n
    
    s1 = f(t(i), w(i));
    s2 = f(t(i)+h2, w(i)+h2*s1);
    s3 = f(t(i)+h2, w(i)+h2*s2);
    s4 = f(t(i)+h, w(i)+h*s3);
    
    w(i+1) = w(i) + h/6*(s1 + 2*s2 + 2*s3 + s4);
end

end


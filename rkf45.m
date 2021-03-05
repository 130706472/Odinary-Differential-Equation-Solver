%% Runge_Kutta_Fehlberg Method Function
function [t,y] = rkf45(fun, tspan, y0, h, rTol)
%{
Spring 2020
AERO 300
Lab7
Eddie Hsieh
Runge_Kutta_Fehlberg method function
------------------------------------
Input:
    m: number of equations
    fun: function to be solved
    tspan: time interval, [t1 t2]
    y0: initial value vector
    h: initial step size
    rTol: error tolerance
Output:
    t: time vector
    y: solutions for ODE
%}

% construct time interval
t(1) = tspan(1);
T = tspan(2);
% H step size and vector to store initial step size
H = [];
H = h;
% first term of w(i)
m = size(y0,2);
w = zeros(1,m);
w(1,:) = y0;
% first term of error
err(1,1:m) = 0.1;
% start iteration i = 1
i = 1;


while T - t(i) > 0
%   Start with initial step size
    h = H(1);
%   s_values
    s1 = [fun(t(i),w(i,:))]';
    s2 = [fun(t(i) + 0.25*h, w(i,:) + 0.25*h*s1)]';
    s3 = [fun(t(i) + (3/8)*h, w(i,:) + (3/32)*h*s1 + (9/32)*h*s2)]';
    s4 = [fun(t(i) + (12/13)*h, w(i,:) + (1932/2197)*h*s1 ...
        - (7200/2197)*h*s2 + (7296/2197)*h*s3)]';
    s5 = [fun(t(i) + h, w(i,:) + (439/216)*h*s1 ...
        - 8*h*s2 + (3680/513)*h*s3 - (845/4104)*h*s4)]';
    s6 = [fun(t(i) + 0.5*h, w(i,:) - (8/27)*h*s1 ...
        + 2*h*s2 - (3544/2565)*h*s3 + (1859/4104)*h*s4 - (11/40)*h*s5)]';
    
%   wtrail and ztrail  
    w(i+1,:) = w(i,:) + h*((25/216)*s1 + (1408/2565)*s3 + ...
        (2197/4104)*s4 - (1/5)*s5);
    z(i+1,:) = w(i,:) + h*((16/135)*s1 + (6656/12825)*s3...
        + (28561/56430)*s4 - (9/50)*s5 + (2/55)*s6);
    
%   relative error for step size control
    err(i+1,:) = abs(z(i+1,:) - w(i+1,:));
%   relative tolerance
    RE(i,:) = err(i+1,:).*(abs(w(i+1,:)).^-1);
    
%   tell wether the relative tolerance larger than rTOL or not   
    if RE(i,:) > rTol         
%       re-size the step size
        hs = 0.8*(((rTol*abs(w(i+1,:)).*(err(i,:).^-1))).^(1/5))*h;
%       pick the smallest step from three reduced stepsize and store to H
        H(end+1) = min(hs);
%       let step size be the latest stored step size
        h = H(end);
        
%       re-apply scheme
%       s_values
        s1 = [fun(t(i),w(i,:))]';
        s2 = [fun(t(i) + 0.25*h, w(i,:) + 0.25*h*s1)]';
        s3 = [fun(t(i) + (3/8)*h, w(i,:) + (3/32)*h*s1 + (9/32)*h*s2)]';
        s4 = [fun(t(i) + (12/13)*h, w(i,:) + (1932/2197)*h*s1 ...
        - (7200/2197)*h*s2 + (7296/2197)*h*s3)]';
        s5 = [fun(t(i) + h, w(i,:) + (439/216)*h*s1 ...
            - 8*h*s2 + (3680/513)*h*s3 - (845/4104)*h*s4)]';
        s6 = [fun(t(i) + 0.5*h, w(i,:)...
            - (8/27)*h*s1 + 2*h*s2 - (3544/2565)*h*s3 + ...
            (1859/4104)*h*s4 - (11/40)*h*s5)]';
    
%       wtrail and ztrail  
        w(i+1,:) = w(i,:) + h*((25/216)*s1 + (1408/2565)*s3 +...
            (2197/4104)*s4 - (1/5)*s5);
        z(i+1,:) = w(i,:) + h*((16/135)*s1 + (6656/12825)*s3 +...
            (28561/56430)*s4 - (9/50)*s5 + (2/55)*s6);

%       calculate error for step size control and relative tolerance again
        err(i+1,:) = abs(z(i+1,:) - w(i+1,:));
        RE(i,:) = err(i+1,:).*(abs(w(i+1,:)).^-1);

%       Tell if the relative tolerance passed or not           
            while RE(i,:) > rTol 

%               reduce step size and apply scheme until pass the test
                h = h/2;
%               s_values
                s1 = [fun(t(i),w(i,:))]';
                s2 = [fun(t(i) + 0.25*h, w(i,:) + 0.25*h*s1)]';
                s3 = [fun(t(i) + (3/8)*h, w(i,:) + (3/32)*h*s1 +...
                    (9/32)*h*s2)]';
                s4 = [fun(t(i) + (12/13)*h, w(i,:) + (1932/2197)*h*s1 ...
                - (7200/2197)*h*s2 + (7296/2197)*h*s3)]';
                s5 = [fun(t(i) + h, w(i,:) + (439/216)*h*s1 ...
                    - 8*h*s2 + (3680/513)*h*s3 - (845/4104)*h*s4)]';
                s6 = [fun(t(i) + 0.5*h, w(i,:) - (8/27)*h*s1 ...
                    + 2*h*s2 - (3544/2565)*h*s3 + (1859/4104)*h*s4 ...
                    - (11/40)*h*s5)]';

%               wtrail and ztrail  
                w(i+1,:) = w(i,:) + h*((25/216)*s1 + (1408/2565)*s3...
                    + (2197/4104)*s4 - (1/5)*s5);
                z(i+1,:) = w(i,:) + h*((16/135)*s1 + (6656/12825)*s3 +...
                    (28561/56430)*s4 - (9/50)*s5 + (2/55)*s6);

%  calculate error for step size control and relative tolerance again
                err(i+1,:) = abs(z(i+1,:) - w(i+1,:));
                RE(i,:) = err(i+1,:).*(abs(w(i+1,:)).^-1);                
            end
    end
    
% let w = z after passed the test    
w(i+1) = z(i+1);
% h is now redused so t(i+1) = t(i) + h_bar, h_bar  = H(end)
t(i+1) = t(i) + h;
% iteration + 1       
i = i + 1;
end
% adjst the last term to step exactly in the time interval
t(i) = T;
% last step size is the end time - previous time, to step on the end
% exactly
h = T - t(i-1);
% extract value from one iteration before the last one
i = i - 1;
% apply scheme
s1 = [fun(t(i),w(i,:))]';
s2 = [fun(t(i) + 0.25*h, w(i,:) + 0.25*h*s1)]';
s3 = [fun(t(i) + (3/8)*h, w(i,:) + (3/32)*h*s1 + (9/32)*h*s2)]';
s4 = [fun(t(i) + (12/13)*h, w(i,:) + (1932/2197)*h*s1 ...
    - (7200/2197)*h*s2 + (7296/2197)*h*s3)]';
s5 = [fun(t(i) + h, w(i,:) + (439/216)*h*s1 ...
    - 8*h*s2 + (3680/513)*h*s3 - (845/4104)*h*s4)]';
s6 = [fun(t(i) + 0.5*h, w(i,:) - (8/27)*h*s1 ...
    + 2*h*s2 - (3544/2565)*h*s3 + (1859/4104)*h*s4 - (11/40)*h*s5)]';

% w trail and z trail  
w(i+1,:) = w(i,:) + h*((25/216)*s1 + (1408/2565)*s3 + (2197/4104)*s4...
    - (1/5)*s5);
z(i+1,:) = w(i,:) + h*((16/135)*s1 + (6656/12825)*s3 + (28561/56430)*s4...
    - (9/50)*s5 + (2/55)*s6);
% last term of w equal to extrapolated term
w(end,:) = z(end,:);

% transpose t to column vector
t = t';
% assign w to y
y = w;
end


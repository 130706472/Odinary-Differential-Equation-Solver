function [t, w] = RK2_3(f, t0, y0, T, TOL, h)

% Implements the 2/3 Embedded Runge-Kutta pair (for ETM and Simpson's method) for 1st-order IVPs
%   f: right-hand side of ODE
%   t0: Initial time
%   y0: Function value at the initial time
%   T: final time at which estimate is desired
%   TOL: prescribed relative tolerance
%   h: prescribed initial step size

err = []; %vector to store LTE at each time step
err(1,1) = 0.1; %Small number instead of true value of the LTE which is zero at starting point

H = []; %vector to store final time step after succesful iteration
H(1,1) = h; %initialize using prescribed initial time step

t(1,1) = t0; %set initial time
w(1,1) = y0; %set value at initial time

while T - t(end,1) > 0.0
    
    % Apply scheme
    
    s1 = f(t(end,1), w(end,1));
    s2 = f(t(end,1) + h, w(end,1) + h*s1);
    s3 = f(t(end,1) + h/2, w(end,1) + h/2*(s1 + s2)/2);
    
    w_trial = w(end,1) + h/2*(s1 + s2);
    z_trial = w(end,1) + h/6*(s1 + 4*s3 + s2);
    
    err_trial = abs(w_trial - z_trial);
    relErr_trial = err_trial/abs(w_trial);
    
    % Enter if error criterion not met
    
    if relErr_trial > TOL
        
        % Compute "optimal" h
        
        h = 0.8*H(end,1)*(TOL*abs(w_trial)/err(end,1))^(1/3);
        
        % Re-apply scheme, skipping s1 as it is independent of h
        
        s2 = f(t(end,1) + h, w(end,1) + h*s1);
        s3 = f(t(end,1) + h/2, w(end,1) + h/2*(s1 + s2)/2);
        
        w_trial = w(end,1) + h/2*(s1 + s2);
        z_trial = w(end,1) + h/6*(s1 + 4*s3 + s2);
        
        err_trial = abs(w_trial - z_trial);
        relErr_trial = err_trial/abs(w_trial);
        
        % Enter loop if using "optimal" h did not work
        
        while relErr_trial > TOL
            
            h = h/2; %repeatedly divide by 2 until convergence
            
            s2 = f(t(end,1) + h, w(end,1) + h*s1);
            s3 = f(t(end,1) + h/2, w(end,1) + h/2*(s1 + s2)/2);
        
        w_trial = w(end,1) + h/2*(s1 + s2);
        z_trial = w(end,1) + h/6*(s1 + 4*s3 + s2);
            
            err_trial = abs(w_trial - z_trial);
            relErr_trial = err_trial/abs(w_trial);
        end
        
    end
    
    w(end+1,1) = z_trial; %extrapolate by using approximation from z
    t(end+1,1) = t(end,1) + h;
    
    H(end+1, 1) = h;
    err(end+1,1) = err_trial;
    
    h = 2.0*H(end, 1); %perhaps increase h so as not to use smaller h if not needed?
end

%As written, algorithm overshoots right end point.
%The following corrects for that.
%Error criterion below is automatically met because h is less than previously computed value

h = T - t(end-1, 1);

s1 = f(t(end-1,1), w(end-1,1));
s2 = f(t(end-1,1) + h, w(end-1,1) + h*s1);
s3 = f(t(end-1,1) + h/2, w(end-1,1) + h/2*(s1 + s2)/2);

z_trial = w(end-1,1) + h/6*(s1 + 4*s3 + s2);

w(end,1) = z_trial;
t(end,1) = T;

end


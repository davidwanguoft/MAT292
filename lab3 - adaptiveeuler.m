function AE_out = adaptiveEuler(DE, t0,tN,y0,h) 

% Initialize arrays
ftemp(1) = y0;
y(1) = y0;
t_array(1) = t0;

% Initialize values
t = t0;
firststepsize = h;
previous = 0; % Counter: revious dy/dt initially 0
n = 1;

% Iterate IEM
while (t + h <= tn) 
    
    % This is the next value of the function 
    ftemp = [ftemp, ftemp(n) + DE(t(n), ftemp(n))*h]; 
    
    current = ((f(t + h, y(n)) - f(t, y(n))) + (f(t, ftemp(n+1)) - f(t, y(n)))* f(t, y(n))) / h;
    
    if (((current > (1.1*previous))) || first_derivative == 0)
        h = h/2; % Step size 1/2 of original
        current = ((f(t + h,y(n)) - f(t,y(n))) + (f(t,ftemp(n+1)) - f(t,y(n)))* f(t,y(n))) / h;
        previous = current;
        y = [y, y(n) + DE(t(n), y(n))*h];
        tempy = y;
    end
    
    y = tempy;
    
    %Next t value
    t_array(n+1) = t_array(n) + h;
    
    h = firststepsize; 
    t = t + h;
    n = n + 1;
end

AE_out = [t_array; y]; %Function of t and y into IE_out

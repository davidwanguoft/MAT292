function improved_euler = IE(DE, t0, tn, y0, h) 

numofsteps = (tn - t0)/h; 

%Initializing the arrays
y = zeros(1, numofsteps + 1); % +1 for of points at each side of step (e.g. split rectangle into 4 pieces with 3 lines)
y(1) = y0;
t = zeros(1, numofsteps + 1);
t(1) = t0;

for n = 1: stepno 
    
    t(n+1) = t(n) + h;                      % Define new t value        
    y(n+1) = y(n) + DE(t(n), y(n))*h;       % Find first y value, nth one
    
    y(n+1) = y(n) + [( DE(t(n), y(n)) + DE(t(n+1), y(n+1))) / 2] * h;       % Use y(n+1) to get IMPROVED Euler.   
end
IE_out = [y; t];
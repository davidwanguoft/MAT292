%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Heun method (or 
% Improved Euler method), and compare its results to those of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in at the
% end of the lab.  Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online using Blackboard.
%
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: David Wang
%
% Student Number: 1001XXXXXX
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now matlab can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check MATLAB lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

% function improved_euler = IE(DE, t0, tn, y0, h) 
%
% numofsteps = (tn - t0)/h; 
%
% Initialize arrays
%
% y = zeros(1, numofsteps + 1); % +1 for of points at each side of step (e.g. split rectangle into 4 pieces with 3 lines)
% y(1) = y0;
%
% t = zeros(1, numofsteps + 1);
% t(1) = t0;
% 
% for n = 1: stepno 
%    
%     t(n+1) = t(n) + h;                      % Define new t         
%     y(n+1) = y(n) + DE(t(n), y(n))*h;       % Find first y, nth one
%    
%     y(n+1) = y(n) + [( DE(t(n), y(n)) + DE(t(n+1), y(n+1))) / 2] * h;       % Use y(n+1) to get IMPROVED Euler.   
% end
%
% IE_out = [y; t];

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the initial-value problems of exercises 1, 3, 4 and 5 
% from MATLAB lab 2, approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the grahs of your Improved Euler Approximation with the |ODE45| 
% approximation.
%
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences of
% note exercise by exercise.


% - ode45 chooses step sizes for you
% - improved euler method makes you choose the step sizes
% - the largest error in ode45 < [less than] that of the euler function (step size matched with number of points)

%% Exercise 3
%
% Objective: Use Euler method and verify that the estimate for the global
% truncation error studied in class is valid.
%
% Details: 
%
% (a) Use Euler method (you can use
% euler.m from iode) to solve the IVP
%
%  |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.75|.
%
% (b) Calculate the solution of the IVP and evaluate it at |t=0.75|.
%
% (c) In lecture, we learned how to estimate the global truncation 
% error for the Euler method. Write the estimate (as a comment), specifying
% the meaning of all the constants used.
%
% (d) Compute the error estimate for |t=0.75| and compare with the actual
% error.
%
% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the global truncation error
% studied in lectures.

% (a)

% f = @(t, y) 2*t*(1 - y^2)^(0.5);
% inline it (http://www.mathworks.com/help/matlab/ref/inline.html)

euler(inline('2.*x.*(1-y^2).^(0.5)','x','y'), 0, linspace(0, 0.75, 100));

% (b)
y(t) = sin(t^2);
y(0.75)

% (c)

% We define error as:
% error <= ((K * h) / (2 * M)) * (exp(M * h * n) - 1),

% h = step size
% M = bound of the partial of f w.r.t t
% f = function (w.r.t y)
% K = bound of double derivative of: K = M + M^2


syms t;
syms y;
f = @(t,y) 2*t*(1 - y^2).^0.5;

fminbnd(f, 0, 0.75);

% (d) 

% globalerror = ((k*timestep)/(2*m)*(exp(m*(tn-t0)) - 1)

% Intialize
syms t tex y

% Set the time step
timestep = 0.01

% Function
f = @(t, y) 2*t*(1 - y^2)^(0.5);
yeu = euler(f, 0, 0: timestep: 0.75)


% Exact solution
exact = sin(tex^2);
exactsoln = eval(subs(exact, tex, 0:timestep:0.75));

yexeval = eval(subs(exact, tex, 0:timestep:0.75));

% Partial w.r.t t
f_t = eval(subs(diff(f, t), t, 0:timestep:0.75));
f_tvals = eval(subs(f_t, yexeval));

% Partial w.r.t y
f_y = eval(subs(diff(f, y), t, 0:timestep:0.75)); 
f_yvals = eval(subs(f_y, yexeval));

% Subbed
ts = eval(subs(f, t, 0:timestep:0.75))
ys = eval(subs, ts, y, yexeval); 

% Compute values
k = max([abs(f_tvals), abs(f_yvals), abs(ys)]);
m = max([abs(f_yvals)]);

esterror = ((k*timestep)/(2*m))*(exp(m*(0.75) - 1))
accerror = abs(exactsoln - y);
actualerror = errAc(end);

% X ((0.946 * (0.5)) / (2 * 2)) * (exp(2 * (0.5) * (0.75/(0.5))) - 1)
% redefine
% rat 1-1
% 

% (e)

% Estimated error: 0.0693
% Actual error: 0.0285
% Confirmed global trunc. error of Euler = 1
% If step size is X times larger, then error is X times larger. 
% I.e. 5-5/4-4/2-2

%% Adaptive Step Size
%
% As mentioned in MATLAB assignment 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where h is the maximum step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Use Euler's method from |iode| and change it to include the following:
%
% (a) Approximation of the derivative of |f(t,y(t))|:
%
%  |d/dt ( f(t,y(t) ) = [ df/dt + df/dy . f ] (t,y(t))|
%
% by 
%
%  |{ [ f(t[n+1],y[n]) - f(t[n],y[n]) ] + [ f(t[n],y[n+1]) - f(t[n],y[n]) ] . f(t[n],y[n]) } / h|
%
% (b) At every iteration compute this approximation and check whether or
% not it is greater than 1.1 times that of the previous iteration.
%
% (c) If it is, then re-compute that iteration with |h/2| instead of |h|
% (decrease the time step by half).
%
% (d) Continue until the approximation reaches the predetermined |tN|.

function[tae, solution] = AE(f, t0, tn, y0, h)

% Assignments
% y = y0;
% solution = [y0];
% tae = [t0];
% i = 1;
% DE = zeroes(size(tae));
% 
% while tae[i] < t                
%     
%     t = tae(i) + h;                             
%     k = feval(f, tae[i], y);                    
%     
%     fex= f(t,soln(i));                          % Soln actual
%     fae = f(tae(i),soln(i));                    % Soln AE
%     faeincr = (f(tae(i), soln(i) + h*k);           % Soln AE with increment
%     
%     DE(i) = (ex - fae) + faeincr - fae) * fae) / h;      % Derivative function (approximate)
%     
%     if (i==1)
%         tae = [tae, t];                         % 
%         y = y + h*k;                            % Reassignment for case 1
%     	  soln=[soln,y];                          % 
%     else
%         while (DE(i) >= (1.1) * DE(i - 1) && h > (10^-5))         % Condition (b)
%             h = 0.5 * h;                          % Half
%             t = tae(i) + h;                     % Step up 
%             
%             DE(i) = (ex - fae) + faeincr - fae)* fae) / h; % Adapt to new initialized values
%                                                            % SAME AS FUNCTION DE ABOVE. REINITIALIZE  
%         end
%         
%         y = y + h*k;                            % |
%         tae = [tae, t];                         % | Reassignment for iteration
%         solution = [solution, y];               % |
%         
%     i = i + 1;                                  % Re-run (iterate)
% 
%     end
  
        
            
%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
%
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
%
% (c) Plot both approximations together with the exact solution.

% (a)

x = linspace(0, 0.75, 2);
y = euler(inline('2.*x.*(1-y^2).^(0.5)','x','y'), 0, linspace(0, 0.75, 2));

% (b)

f = @(t,y) 2*t*(1 - y^2)^(0.5);
solution = adaptiveeuler(f, [0, 0.75], 0, 0.025);

% (c)

figure;
plot(solution.x, solution.y, '+', 'MarkerSize', 10, 'MarkerEdgeColor', 'red', 'LineWidth', 2);
hold on;

plot(x, y, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
grid on;

hold off; 


%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from Exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% 
% (b) Plot the Exact Solution (done in 3.b), the Euler's Approximation
% (done in 3.a) and the Adaptive Euler's Approximation (done in 5) from
% |t=0| to |t=1.5|.
%
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

clear all;

% (a) 
% The approximations are the same, for observations derrived from the
% graph. The slope decreases over t++, and will never be greater than 1.1
% times the previous slope value. 
% The defined condition means they will be identical.

% (b)

% Euler
f = @(t,y) 2*t*(1 - y^2)^(0.5);

ye = euler(inline('2.*x.*(1-y^2).^(0.5)','x','y'), 0, linspace(0, 0.75, 100))

% Adaptive Euler

[ta, yae] = adaptiveeuler('2.*x.*(1-y^2).^(0.5)', 0, 1.5, 0, 0.025);

% Exact

yex = sin(t.^2);

plot(t, ye, ta, yae, t, yex);
xlabel('t');
ylabel('y');
legend('Euler', 'Adapt Euler', 'Exact', 'Location', 'Best');

% (c)

% The exact solution and appeoximations become very different and an
% error ensues. There is a square root in the function, and for the
% function to be defined, C in sqrt(C) must be real. If y > 1, then the
% root is imaginary.

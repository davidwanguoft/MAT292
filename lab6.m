%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|. Also in this lab, you will write your own
% ODE solver using Laplace transforms and check whether the result yields
% the correct answer.
%
% You will learn how to use the |laplace| routine. 
% 
% There are five (5) exercises in this lab that are to be handed in at the
% end of the lab.  Write your solutions in the template, including
% appropriate descriptions in each step. Save the m-file and submit it 
% online using Blackboard.
%
% Include your name and student number in the submitted file.
%
% MAT292, Fall 2015, Sousa, based on 
% MAT292, Fall 2013, Sinnamon & Sousa

%% Student Information
%
%  Student Name: David Wang
%
%  Student Number: 1001XXXXXX
%

%% Using symbolic variables to define functions
% 
% Recall the use of symbolic variables and function explained in the MATLAB
% assignment #2.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
% 
% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.

syms t s

% Part A
% Defining function f(t) and computing it's laplace transform
f = exp(2*t)*t^3;
F = laplace(f);
display(F);

% Part B
% Defining function G(s), which the laplace transform of a function g(t).
q = ((s - 1)*(s - 2))/(s*(s + 2)*(s - 3));
% Calculating the inverse laplace transform of G(s) to find g(t)
Q = ilaplace(q);
display(Q);

% Part C
% Defining a function g(t) which does not have a laplace transform
g = 1/sqrt(t^2+1);
% Trying to compute the laplace transform of g(t) using |laplace| so MATLAB
% knows that the variable G is the laplace transform of g(t)
G = laplace(g);
% Calculating the laplace transform of exp(2*t)*g(t), however we can use
% |ilaplace| to give us g(t). The answer shows that laplace transform is
% computed by finding the laplace transform of g(t) but with the variables
% t and (s-2). 
Gtransform = laplace(exp(2*t)*ilaplace(G));
display(Gtransform); 




%% Heaviside and Dirac functions
%
% These two functions are builtin in MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function is defined at |0| also with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)

%% Caution! 
%
% MATLAB doesn't handle well function multiplying the dirac delta.
%
% Consider a function |f(t)| and let |f2 = f(2)|.

syms f f2

g = dirac(t-3)*f

laplace(g)

% Compare this with the correct result

G = exp(-2*s)*f2


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% (in your answer, explain the 'proof' using comments).

syms t f(t)
a = 2;

g = heaviside(t-a)*f(t-a);
f = f(t);

G = laplace(g);
F = laplace(f);

display(G);
display(F);

% If G(s) = laplace(u_a(t)f(t-a)), and F(s) = laplace(f(t))
% Note that the only difference is that G(s) is the laplace transform
% computed for f(t) translated by t-a
% Then G(s) = exp(-a*s)*F(s)



%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t 

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,'y(0)',1)
L_ODE=subs(L_ODE,'D(y)(0)',2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,'laplace(y(t), t, s)', Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: (Explain your steps using comments)
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)+100*dirac(t-50*pi)|
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
%
% * Plot the solution for |t| in |[0,300]|
%
% * Without the |dirac| function, the system resonates. Comment on how the
% |dirac| affects the resonance.
%

syms y(t) t Y

DE = diff(y(t),t,3) + 2*diff(y(t),t,2) + diff(y(t),t,1) + 2*y(t) + cos(t) - 100*dirac(t-50*pi);
L_DE = laplace(DE);

L_DE = subs(L_DE, 'y(0)', 0);
L_DE = subs(L_DE, 'D(y)(0)', 0);
L_DE = subs(L_DE, 'D(D(y))(0)', 0);

L_DE = subs(L_DE, 'laplace(y(t), t, s)', Y);
lTransform = solve(L_DE, Y);

y = ilaplace(lTransform);

ezplot(y, [0,300]);





%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
%  (in your answer, explain your steps using comments).

syms y(t) t Y 
g = 3*heaviside(t)+heaviside(t-2)*(t+1-3)+heaviside(t-5)*(5-t-1);
fDE = diff(y(t),t,2) + 2*diff(y(t),t) + 5*y(t) - g;

LTDE = laplace(fDE);

LTDE = subs(LTDE, 'y(0)', 2);
LTDE = subs(LTDE, 'D(y)(0)', 1);
LTDE = subs(LTDE, 'laplace(y(t), t, s)', Y);

LTDE = solve(LTDE, Y);

y = ilaplace(LTDE);

q = ezplot(y, [0,12, 0, 2.25]);

% Uncomment if you're using MatLab 2013
%X = get(q, 'XData');
%Y = get(q, 'YData');

%plot(X,Y)
%set(gca, 'XLim', [0,12], 'YLim', [0,2.5]);
% Uncomment if you're using MatLab 2013


%% Exercise 5
%
% Objective: Solve an IVP with a periodic function using the Laplace
% transform.
%
% Details:
%
% * Follow the instructions on the website
%  http://instruct.math.lsa.umich.edu/lecturedemos/ma216/docs/7_5/
%
% * Do the part labelled |Outside of Lecture|
%
% * Note that |u(t-a)| is the Heaviside function |u_a(t)| that we
% defined in lecture
%
% * Check Theorem 5.5.3 (page 349) to know how to define the Laplace 
% transform of a periodic function like the one in this exercise (and check
% the function |int| on MATLAB for symbolic integration).
%
% * You are using the only second-order DE there
% * You should use the numbers provided (the less precise ones)
% * You should obtain an exact solution. In the process, you should need symbolic integrations.
%
% * The goal is to solve the problem described (it also shows what your
% solution should look like).


syms t f Y s y(t) g

% Calculate LaPlace of Forcing Function
f = 0.48*(heaviside(t)-heaviside(t-6));
f = int(exp(-s*t)*f,t,0,24);
f = f/(1-exp(-s*24));

DE = diff(y(t),t,2) + 0.54*diff(y(t),t,1) + 0.02*y(t);

L_DE = laplace(DE) - f;

L_DE = subs(L_DE,'y(0)',0);
L_DE = subs(L_DE,'D(y)(0)',0);
L_DE = subs(L_DE,'laplace(y(t), t, s)', Y);
L_DE = solve(L_DE, Y);

g = ilaplace(L_DE);

ezplot(g,[0, 250, 0, 10]);




%------------------------ARC-LENGTH CONTINUATION METHOD--------------------
%The arc-length continuation method is a numerical technique used to trace 
% the response of nonlinear equations while maintaining stability and 
% accuracy. It is particularly useful in situations where the 
% Newton-Raphson method encounters difficulties, such as when the equations
% exhibit snap-back behavior in their plots.

%In the arc-length continuation method, an additional equation in the form 
% of a circle equation is introduced to the original set of nonlinear 
% equations. This equation provides an arc-length constraint, which allows 
% for controlled parameter variation along the response curve. 
% By incorporating the arc-length equation, the path of the nonlinear 
% equations can be followed smoothly and accurately.

%To apply the arc-length continuation method, the parameter of interest,
% denoted as w, is systematically varied, and the resulting system of 
% equations is solved using the Newton-Raphson method. At each step, 
% the equations are solved iteratively to obtain the corresponding solution
% for the unknown variables, with the arc-length equation serving as an 
% additional constraint. By successively adjusting w and solving the 
% equations, the entire response curve of the nonlinear system can be traced.

%References:

%[1] de Borst, R., Crisfield, M. A., Remmers, J. J. C. & Verhoosel, 
%    C. V. NON-LINEAR FINITE ELEMENT ANALYSIS OF SOLIDS AND STRUCTURES WILEY 
%    SERIES IN COMPUTATIONAL MECHANICS Introduction to Finite Strain Theory 
%    Hashiguchi and Yamakawa September 2012 for Continuum Elasto-Plasticity. 
%    感染症誌 vol. 91 (Wiley, 2012).

%[2] OpenAI. (2021). GPT-3.5 "ChatGPT" [Software]. Retrieved from https://openai.com

%[3] Nayfeh, A. H. & Mook, D. T. Nonlinear Oscillations. Classical Dynamics 
% of Particles and Systems (Wiley, 1995). doi:10.1002/9783527617586.

clc;
clear;
close all

%% nonlinear equations[3]:
alpha = 2;            % nonlinearity coeff
A = 1;                % force amplitude
w0 = 1;               % natural frequency
mio = 0.7;            % damping coeff
syms w                % excitation frequency
n = 1;                % number of unknown variables (degrees of freedom or number of equations)
syms('u',[1 n])       % amplitudes
syms('ustar1',[1 n])  % center of the arc x direction(previous point)
syms R1 beta1 wstar1  % arc radious, user specified parameter, center of the arc y direction

%equations:
%insert equations f = [f1(u1,u2,...,w);f2(u1,u2,...,w);...]
f = [(mio^2 + (w - (3*alpha*A^2)/w0 - (3*alpha*u1^2)/(8*w0))^2)*u1^2 - (alpha^2*A^6)/(w0^2)];
%% define curve equation:
u = u.';
ustar1 = ustar1.';

%equation of the curve:
curveq = ([u]-ustar1).'*([u]-ustar1)+beta1*(w-wstar1)^2-R1^2;

eqns = [f;curveq];

%% solve using fsolve
% choose the first point:
w1 = 0.1;       %first point
w2 = 0.2;       %second point

Uin = zeros(n,1); %initial guess
%solving the nonlinear equations for two small number of frequency in order
%to find the first point on the curve and also the direction of the curve:
nonlinear_function1 = matlabFunction(subs(f,w,w1),'Vars',{u});
nonlinear_function2 = matlabFunction(subs(f,w,w2),'Vars',{u});
U01 = fsolve(nonlinear_function1,Uin);
U02 = fsolve(nonlinear_function2,Uin);

%ARC LENGTH CONTINUATION:
beta = 0.001;                %a user-specified value that weighs the importance of the contributions
u10 = [U02;w2];            %first point on the curve
dir = ([U02;w2]-[U01;w1]); %direction of the curve

ustar = U02;               %center of the arc
wstar = w2;                %center of the arc
ustar0 = U02;
wstar0 = w2;
R = 0.01;          %arc radius

figure;
ax = gca;
lineObj = animatedline(ax);

ustar = [U02;w2];
for nit = 1:2000 %numbear of iteration that ARC-LENGTH CONTINUATION must follow the path


arclength = matlabFunction(subs(eqns,[ustar1.',wstar1,beta1,R1],[ustar.',beta,R]),'Vars',{[u;w]});
disp("number of iteration = " + num2str(nit))

    if nit == 2
        dir = (u_sol-[ustar0;wstar0]);     %direction of the curve
        u10 = u_sol(:,nit-1);              %center of the arc
    elseif nit>2
        u10 = u_sol(:,nit-1);              %center of the arc
        dir = (u_sol(:,nit-1)-u_sol(:,nit-2));%direction of the curve
    end

    u_sol(:,nit) = fsolve(arclength,[u10+(dir)]);

    drawnow;
amplitude(:,nit) = (u_sol(1:end-1,nit));
frequency(nit) = double(u_sol(end,nit));
required_sol = amplitude(1,nit);
addpoints(lineObj,frequency(nit),required_sol)
drawnow;

% Adjust the axis limits to zoom in on the new point
xlim([frequency(nit) - 0.1, frequency(nit) + 0.1]);
ylim([required_sol - 0.1, required_sol + 0.1]);

ustar = u_sol(:,nit);
wstar = u_sol(end,nit);
end


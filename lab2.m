%% [Part 1] Simulate the dynamics of the passive membrane model with different Vm(0)

%Declare constants
C_m = 3; % Membrane capacitance = 3 
G_L = 4; % Leaky conductance = 4
E_L = 1; % Resting membrane potential = 1mV

%time conditions and steps (range from 0~10 seconds)
t0 = 0; 
tf = 10; 
dt = 0.001; %step size 1ms

%time vector and voltage Vector
tvec = t0:dt:tf;
vvec = zeros(size(tvec));

%declare initial voltage vector
Vm0_initial = zeros(3); 
Vm0_initial(1) = 2; % Vm0>EL
Vm0_initial(2) = 1; % Vm0=EL
Vm0_initial(3) = 0.5; % Vm0<EL

% Implementing Euler Methods
for j = 1:3
    %implementing different Vm0 initial conditions
    vvec(1) = Vm0_initial(j);
    %Euler Methods
    for i = 2:length(tvec)
    %dvm/dt = Gl/Cm*(E_L-Vm) 
    dVm_over_dt =  G_L / C_m *(E_L - vvec(i-1));
    % Vm(i) = Vm(i-1) + dt*f(t(i-1),Vm(i-1)) 
    vvec(i) = vvec(i-1) + dt * dVm_over_dt;
    end
    figure(1)
    plot(tvec, vvec)
    hold on
end 
%labels, legends and axis
title('membrane potential vs time with different Vm0s')
xlabel('time (ms)')
ylabel('membrane potential (Vm)')
legend({'Vm(0)>E_L','Vm(0)=E_L','Vm(0)<E_L'})

%% What can you say about the dynamics of this model behavior?

% When the initial Vm(0) equal to Equilibrium potential, dVm/dt = 0.
% Therefore, we can see a constant line parallel to the x axis. However, 
% when Vm(0) is greater than EL, dVm/dt turns negative, which can
% be seen a negative exponential graph Vm approaching EL as time goes on.
% On the other hand, when Vm(0) is less than EL, dVm/dt is positive, which
% translates to a postive exponential graph which Vm approaches to EL over
% time.

%% [Part 2] Simulate the coupled two variable system of equations on Slide 9

% define initial constants 
x0 = 1.0; % x(0)=1m
v0 = 0.0; % v(0)=0m/s
omega = 3; % omega = 3
omega_sqr = omega * omega; % omega^2 = 9

%time vector from 0s to 100s
t0 = 0;
tf = 100;
dt = 0.00001; 
tvec = t0:dt:tf;

%define the velocity and position vector
vvec = zeros(size(tvec));
xvec = zeros(size(tvec));

%implement initials from the velocity and position vector
vvec(1) = v0;
xvec(1) = x0;

%Euler's Method implementation.
for i = 2:length(tvec)
   %dx(t)/dt = v(t)
   dx_over_dt = vvec(i-1);
   %dv(t)/dt = -omega^2*x(t)
   dv_over_dt = -1 * omega_sqr * xvec(i-1);
   %yi = yi-1 + dx*f(xi-1,yi-1)
   xvec(i) = xvec(i-1) + dx_over_dt *dt;
   vvec(i) = vvec(i-1) + dv_over_dt *dt;
end

%plot the solutions
%plot 1
figure(2)
plot(tvec, xvec)
hold on
plot(tvec, vvec)
title('Part 2: oscillation fuctions')
xlabel('time(t)')
ylabel('values (x(t),v(t))')
legend({'dx/dt = v(t)','dv/dt = -omega^2*x(t)'})
%plot 2
figure(3)
plot(xvec,vvec)
title('Part 2: v(t) vs. x(t)')
xlabel('x(t)')
ylabel('v(t)')

%% What can you say about the model behavior?
% we can see that both graphs are sinusoidal. These two models affect each
% other, as can be seen in the circle in the graph. In polar coordinates
% coupled oscillators represent a sine cosine relationship.










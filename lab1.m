%[30 points] submit MATLAB codes from the following.

%problem 1
q = (100:-1:1);

%problem 2
input = [8,8,8,8,8,8,8,8,8,8];
y = -3 * ones(10,10) + diag(input);
disp(y);

%problem 3
z = randi([10,20],5,3);
disp(z);

%problem 4
a = [rand(1,15) + 2 rand(1,5) - 1];
disp(a);

%problem 5
x2 = linspace(-5,5,100);
y1 = 3*x2;
plot(x2,y1);
hold on
y2 = 4*x2;
plot(x2,y2);

%problem 6
 run problem6.m 
%the file is being put in the same project

% problem 7
tot = 0;
for i = 1:15
    %tot = tot + i^3;
end
disp(tot)

%problem 8
for j = 1:999
    if j + j^2 + j^3 + j^4 == 88740
        disp(j)
        break
    end
end

%euler's method

%problem 1
%initial conditions
t0 = 0;
x0 = 0;
%parameters
v0 = 2;
dt = 0.001;
tf = 5;
% Initialize
t = t0:dt:tf;
x = zeros(size(t));
x(1) = x0;
% approximation
for i = 1:length(t)-1
    x(i+1) = x(i) + v0*(1 + x(i)^2)*dt;
end
plot(t, x)

%problem 2
%initial conditions
t0 = 0;
T0 = 20;
%parameters
Te = 25;
tau = 10;
tf = 30;
dt = 0.01;
% Initialize
t = t0:dt:tf;
T = zeros(size(t));
T(1) = T0;
% Euler's method
for i = 1:length(t)-1
    T(i+1) = T(i) + dt*(Te - T(i))/tau;
end
plot(t, T)


%ode45 PROBLEMS

%problem 1
%differential equation
f = @(t, x) v0*(1 + x^2);
% conditions
t0 = 0;
x0 = 0;
%parameters
v0 = 2;
tf = 10;
% ode function
[t, x] = ode45(f, [t0 tf], x0);
plot(t, x)
xlabel('Time')
ylabel('Position')
title('dx/dt = v0*(1 + x(t)^2)')

%problem 2
%differential equation
T1 = 15;
f = @(t, T) (T1 - T)/tau;
% initial conditions
t0 = 0;
T0 = 10;
% parameters
tau = 10;
tf = 10;
% ode45
[t, T] = ode45(f, [t0 tf], T0);
plot(t, T)












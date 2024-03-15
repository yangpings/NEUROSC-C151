%% Reference [Part D] Hodgkin and Huxley Model
function dy = lab5partd(t,y,Iapp) 
Vm = y(1); %set matrix column to Vm
n = y(2); %set matrix column to n
m = y(3); %set matrix column to m
h = y(4); %set matrix column to h

% Setting up Parameters in SI units
Gl = 30*10^-9; %Leak Conductance 30nS
GNa = 12*10^-6; %Max Sodium Conductance 12uS
GK = 3.6*10^-6; %Max Delayed Rectifier 3.6uS
ENa = 45*10^-3; %Sodium Reversal Potential 45mV
EK = -82*10^-3; %Potassium Reversal Potential -82mV
EL = -60*10^-3; %Leak Reversal Potential -60mV
Cm = 100*10^-12;%Membrane Capacitance 100pF
delay = Iapp(1); %setting up delay
%time/applied current
start =zeros(1,10); %applied current start time vector
finish =zeros(1,10); % applied current end time vector
start(1) = 0; %first start time
finish(1) = 0.005; % first end time
for i = 2:10 %loop the values into the vectors
    start(i) = start(i-1) + delay;
    finish(i) = finish(i-1) + delay;
end
for j = 1:10
    if t>=start(j) && t<= finish(j) % if t is within this range
    iapp = 0; %current = 0nA
    break
    else 
    iapp = Iapp(2); % or else it's always 0.6nA
    end
end

%Gating Variable m
alpha_m = (10^5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
beta_m = 4*10^3*exp((-Vm-0.070)/0.018);
%Gating Variable n
if Vm == -0.06
    alpha_n = 100;
else
    alpha_n = (10^4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
end
beta_n = 125*exp((-Vm-0.070)/0.08);
%Gating Variable h
alpha_h = 70*exp(50*(-Vm-0.07));
beta_h = 10^3/(1+exp(100*(-Vm-0.040)));

%joint ODEs
Leaky = Gl*(EL-Vm);
INa = GNa*m^3*h*(ENa-Vm);
IK = GK*n^4*(EK-Vm);
dVm_over_dt = (1/Cm)*(Leaky+INa+IK+iapp);
dn_over_dt = alpha_n *(1-n) - beta_n*n;
dm_over_dt = alpha_m *(1-m) - beta_m*m;
dh_over_dt = alpha_h *(1-h) - beta_h*h;
dy = [dVm_over_dt;dn_over_dt;dm_over_dt;dh_over_dt];
end

%ignore the comments as the function is perfectly integrated into the main
%code without any error message. The function cannot work by its own. 
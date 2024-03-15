
function dy = Morris_Lecar_1(t,y,I)
%HOPF
%Iapp
Iapp_Hopf = I; % in nA

%output Variables (Anthony)
Vhopf=y(1);
nhopf=y(2);

%define Global Variables (Anthony)
Cm = 20 ; %microfarad/cm^2 
ECa=120; %millivolts
gK=8;% millisiemens/ cm^2 
EK=-84; %millivolts
gL=2;% millisiemens/ cm^2 
EL=-60;%millivolts
v1=-1.2; %millivolts
v2= 18 ; %millivolts

%____________________Hopf____________________
%Additional Parameters for Hopf Bifurcation Class 2 Neurons (Anthony)
v3_hopf= 2 ; %millivolts
v4_hopf= 30; %millivoltshb 
phi_hopf=0.04; % per millisecond
gCa_hopf=4.4; % millisiemens/ cm^2 

%Other supporting functions in functions of voltage.
m_inf_hopf = (0.5*(1+tanh((Vhopf-v1)/v2))); %(Anthony)
n_inf_hopf = (0.5*(1+tanh((Vhopf-v3_hopf)/v4_hopf))); %(Anthony)
tau_n_hopf = (1/(cosh((Vhopf-v3_hopf)/(2*v4_hopf)))); % (Eric)

%differential equations for hopf 
% dn/dt (Anthony)
dnhopf_over_dt = phi_hopf*(n_inf_hopf - nhopf) / tau_n_hopf;
% dV/dt (Eric)
P1 = gL*(Vhopf - EL); %leak component
P2 = gK * nhopf * (Vhopf - EK); % pottasium component
P3 = gCa_hopf * m_inf_hopf * (Vhopf-ECa); %Calcium component
dVhopf_over_dt = (Iapp_Hopf - P1 - P2 - P3) / Cm; %combine
%___________________________________________________________

%Wrapper function (Eric)
dy = [dVhopf_over_dt; dnhopf_over_dt];

end
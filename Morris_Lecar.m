
function dy = Morris_Lecar(t,y,I1,I2)
%Iapp
Iapp_Hopf = I1; % in nA
Iapp_SNLC = I2;

%output Variables (Anthony)
Vhopf=y(1);
nhopf=y(2);
VSNLC=y(3);
nSNLC=y(4);

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

%_______________________SNLC functions_______________________
%additional Parameters for (SNLC) Bifurcation Class 1 Neurons (Eric)
v3_SNLC = 12; % in mV
v4_SNLC = 17.4; % in mV
phi_SNLC = 0.067; % in ms^-1
gCa_SNLC=4.4; % millisiemens/ cm^2 

%Other supporting functions in functions of voltage.
m_inf_SNLC = (0.5*(1+tanh((VSNLC-v1)/v2))); %(Anthony)
n_inf_SNLC = (0.5*(1+tanh((VSNLC-v3_SNLC)/v4_SNLC))); %(Anthony)
tau_n_SNLC = (1/(cosh((VSNLC-v3_SNLC)/(2*v4_SNLC)))); % (Eric)

%differential equations for hopf 
% dn/dt(Anthony)
dnSNLC_over_dt = phi_SNLC*(n_inf_SNLC - nSNLC) / tau_n_SNLC;
% dV/dt(Eric)
P1_SNLC = gL*(VSNLC - EL); %leak component
P2_SNLC = gK * nSNLC * (VSNLC - EK); % pottasium component
P3_SNLC = gCa_SNLC * m_inf_SNLC * (VSNLC-ECa); %Calcium component
dVSNLC_over_dt = (Iapp_SNLC - P1_SNLC - P2_SNLC - P3_SNLC) / Cm; %combine
%_____________________________________________________________
%Wrapper function (Eric)
dy = [dVhopf_over_dt; dnhopf_over_dt; dVSNLC_over_dt;dnSNLC_over_dt];

end
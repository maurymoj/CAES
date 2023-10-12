% Steady state analysis of air flow in a pipe for varying lengths
% In the data set lengths vary from 0.330 km up to 140 km and diameters from
% 0.1 to 1.2 m

pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter

%% Initial case - 1 km pipe
% Ambient conditions
T_a = 12+273.15; % T = 12 oC - ~280 K
P_a =101325 ; % P = 101.325 kPa
rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;

% Compressor outlet/pipe inlet
eta_pol = 0.85;

P_1 = 7e6; % 7 MPa
Q_dot = 10; % m3/s
% Q_dot = 34e6/(24*3600) % Conversion from mcmd (millions of cubic meters
% per day to cubic meters per second) - real demand estimation

w_dot = cp*T_a*( (P_1/P_a)^( (gam-1)/(eta_pol*gam) ) - 1);
T_1 = T_a*(P_1/P_a)^( (gam-1)/(eta_pol*gam) )

rho_1 = py.CoolProp.CoolProp.PropsSI('D','P',P_1,'T',T_1,'Air');
nu_1 = py.CoolProp.CoolProp.PropsSI('V','P',P_1,'T',T_1,'Air');
W_dot = rho_1*w_dot;
m_dot = rho_1*Q_dot

% Pipe
D = 0.900; % 900 mm
A = pi*D^2/4;

V_1 = m_dot/(rho_1*A)

Re = rho_1*V_1*D/nu_1








% Actual pipe
L = 1000; % 1 km pipe
eps = 0.04e-3; % Absolute roughness 0.04 mm
epsD = eps/D;

% Iterations for thet estimation of friction factor using Colebrook equation
f = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.
df = 10;
count=0;

while (df > 0.0001 & count < 10) 
    f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2);
    df = (f_n - f)/f;
    f = f_n;
    count = count + 1; 
end
f
% Friction factor for gases
f_W = 0.032*D^-0.333        % Waymouth
f_PanA = 0.085*Re^-0.147    % Panhandle A - Low pressures
f_PanB = 0.015*Re^-0.039    % Panhandle B - High pressures






% P_2 = ?

%% Flow equation "Gas pipeline hydraulics" - m3/day
G = 1; % gas gravity

Qd = 1.1494e-3*(T_a/P_a)*( (P_1^1-P_2^2)/(G*T_f*L*Z*f) )^0.5*D^2.5



%% Sanity check
% Ambient conditions
T_a = 280; % T = 15 oC
P_a = 101325 ; % P = 101.325 kPa

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;

% Compressor outlet/pipe inlet
eta_pol = 0.85;

P_1 = 600000; % 600 kPa
m_dot = 0.02; % kg/s

w_dot = cp*T_a*( (P_1/P_a)^( (gam-1)/ (eta_pol*gam) ) - 1);
T_1 = T_a*(P_1/P_a)^( (gam-1)/ (eta_pol*gam) )

rho_1 = py.CoolProp.CoolProp.PropsSI('D','P',P_1,'T',T_1,'Air');

W_dot = m_dot*w_dot

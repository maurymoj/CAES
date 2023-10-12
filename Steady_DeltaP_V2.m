pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter

%% Initial case - 1 km pipe - Using Panhandle B flow equation

%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15+273.15; % T = 15 oC - ~280 K
P_a = 101325 ; % P = 101.325 kPa
rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
G = 1;      % Specific gas gravity - for air G = 1

P_1 = 7e6; % 7 MPa
T_1 = T_a;

% Flow rate
Q_a = 14*1000000/(24*3600); % Conversion from mcmd (millions of cubic meters
% per day to cubic meters per second) - real demand estimation 70 mscm/day
% St Fergus, 14 is the proportional relative to area of one of the 3, 900
% mm diameter pipes

% Pipe
D = 0.900; % 900 mm
L = 1000; % 1 km pipe
eps = 0.04e-3; % Absolute roughness 0.04 mm

E = 0.8;   % Pipeline efficiency
% Usually varies from 0.6 to 0.92 as a function of liquid presence and age
% "Handbook of natural gas transmission and processing"
%----------------------------------------------------------------------%

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;

rho_1 = py.CoolProp.CoolProp.PropsSI('D','P',P_1,'T',T_1,'Air');
nu_1 = py.CoolProp.CoolProp.PropsSI('V','P',P_1,'T',T_1,'Air');
Z_1 = py.CoolProp.CoolProp.PropsSI('Z','P',P_1,'T',T_1,'Air');

%Q_dot = 10; % m3/s
%m_dot = rho_1*Q_dot;
% Pipe
A = pi*D^2/4;

% Re = rho_1*V_1*D/nu_1
% more adequate expression for Re in pipelines ("Gas pipeline hydraulics")
D_mm = 1000*D;
P_a_kPa = P_a/1000;
Q_a_day = Q_a*24*3600;
nu_1_Po = 10*nu_1;

% Velocity m/s
u_1 = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z_1*T_1/P_1); % Q_b has to be converted to m3/day and D to mm

Re = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_1_Po*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm

% Friction factor
% Iterations for the estimation of friction factor using Colebrook equation
epsD = eps/D;
f = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.
df = 10;
count=0;
while (df > 0.0001 & count < 10) 
    % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
    f_n = (-2*log10(epsD/3.7 + 2.825/(Re*f^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
    df = (f_n - f)/f;
    f = f_n;
    count = count + 1; 
end
f

% Flow equation Q - Standard flow rate m3/day 
% Q = 3.7435e-3*E*(T_a/P_a_kPa)*...
%     ( ( (P_1_kPa^2 - e^s*P_2_kPa^2)/(G*T_f*L_e_km*Z) )^0.5*D_mm^2.667 )  % Weymouth equation
% Q = 4.5965e-3*E*(T_a/P_a_kPa)^1.0788*...
%     ( ( (P_1_kPa^2 - e^s*P_2_kPa^2)/(G^0.8539*T_f*L_e_km*Z) )^0.5394*D_mm^2.6182 )  % Panhandle A equation
% Q = 1.002e-2*E*(T_a/P_a_kPa)^1.02*...
%     ( (P_1_kPa^2 - e^s*P_2_kPa^2)/(G^0.961*T_f*L_e_km*Z) )^0.51*D_mm^2.53 % Panhandle B - full
% L_e and e^s take the elevation change into account
% Q = 1.002e-2*E*(T_a/P_a_kPa)^1.02*...
%     ( (P_1_kPa^2 - P_2_kPa^2)/(T_f*L_km*Z) )^0.51*D_mm^2.53 % Panhandle B - Simplified (No Dz, G=1)
P_1_kPa = P_1/1000;
T_2 = T_1;              % Assuming isothermal flow
T_f = T_2;
P_2_kPa = 0.98*P_1_kPa; % Initial guess
L_km = L/1000;
dP = 10;
count = 0;

while (dP > 0.0001 & count < 10)
    P_f = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
    Z_f = py.CoolProp.CoolProp.PropsSI('Z','P',P_f*1000,'T',T_2,'Air');

    P = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*T_f*L_km*Z_f );
    dP = (P - P_2_kPa)/P_2_kPa;
    P_2_kPa = P;
    count = count + 1;
end

P_2_kPa
P_2 = 1000*P_2_kPa;
Z_2 = py.CoolProp.CoolProp.PropsSI('Z','P',P_2,'T',T_2,'Air');
nu_2_Po = py.CoolProp.CoolProp.PropsSI('viscosity','P',P_2,'T',T_2,'Air');

u_2 = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z_2*T_2/P_2)
Re_2 = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_2_Po*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm

%% Initial case - 1 km pipe using general gas flow equation

%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15+273.15; % T = 15 oC - ~280 K
P_a = 101325 ; % P = 101.325 kPa
rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
G = 1;      % Specific gas gravity - for air G = 1

P_1 = 7e6; % 7 MPa
T_1 = T_a;

% Flow rate
Q_a = 14*1000000/(24*3600); % Conversion from mcmd (millions of cubic meters
% per day to cubic meters per second) - real demand estimation 70 mscm/day
% St Fergus, 14 is the proportional relative to area of one of the 3, 900
% mm diameter pipes

% Pipe
D = 0.900; % 900 mm
L = 70000; % 1 km pipe
eps = 0.04e-3; % Absolute roughness 0.04 mm

%----------------------------------------------------------------------%

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;

rho_1 = py.CoolProp.CoolProp.PropsSI('D','P',P_1,'T',T_1,'Air');
nu_1 = py.CoolProp.CoolProp.PropsSI('V','P',P_1,'T',T_1,'Air');
Z_1 = py.CoolProp.CoolProp.PropsSI('Z','P',P_1,'T',T_1,'Air');

%Q_dot = 10; % m3/s
%m_dot = rho_1*Q_dot;
% Pipe
A = pi*D^2/4;

% Re = rho_1*V_1*D/nu_1
% more adequate expression for Re in pipelines ("Gas pipeline hydraulics")
D_mm = 1000*D;
P_a_kPa = P_a/1000;
Q_a_day = Q_a*24*3600;
nu_1_Po = 10*nu_1;

% Velocity m/s
u_1 = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z_1*T_1/P_1); % Q_b has to be converted to m3/day and D to mm

Re = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_1_Po*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm

% Friction factor
% Iterations for the estimation of friction factor using Colebrook equation
epsD = eps/D;
f = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.
df = 10;
count=0;
while (df > 0.0001 & count < 10) 
    % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
    f_n = (-2*log10(epsD/3.7 + 2.825/(Re*f^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
    df = (f_n - f)/f;
    f = f_n;
    count = count + 1; 
end
f

% Flow equation Q - Standard flow rate m3/day 
% Q = 3.7435e-3*E*(T_a/P_a_kPa)*...
%     ( ( (P_1_kPa^2 - e^s*P_2_kPa^2)/(G*T_f*L_e_km*Z) )^0.5*D_mm^2.667 )  % Weymouth equation
% Q = 4.5965e-3*E*(T_a/P_a_kPa)^1.0788*...
%     ( ( (P_1_kPa^2 - e^s*P_2_kPa^2)/(G^0.8539*T_f*L_e_km*Z) )^0.5394*D_mm^2.6182 )  % Panhandle A equation
% Q = 1.002e-2*E*(T_a/P_a_kPa)^1.02*...
%     ( (P_1_kPa^2 - e^s*P_2_kPa^2)/(G^0.961*T_f*L_e_km*Z) )^0.51*D_mm^2.53 % Panhandle B - full
% L_e and e^s take the elevation change into account
% Q = 1.002e-2*E*(T_a/P_a_kPa)^1.02*...
%     ( (P_1_kPa^2 - P_2_kPa^2)/(T_f*L_km*Z) )^0.51*D_mm^2.53 % Panhandle B - Simplified (No Dz, G=1)
P_1_kPa = P_1/1000;
T_2 = T_1;              % Assuming isothermal flow
T_f = T_2;
P_2_kPa = 0.98*P_1_kPa; % Initial guess
L_km = L/1000;
dP = 10;
count = 0;

while (dP > 0.0001 & count < 10)
    P_f = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
    Z_f = py.CoolProp.CoolProp.PropsSI('Z','P',P_f*1000,'T',T_2,'Air');

    %P = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*T_f*L_km*Z_f );
    P = sqrt(P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*T_f*L_km*Z_f*f );
    dP = (P - P_2_kPa)/P_2_kPa;
    P_2_kPa = P;
    count = count + 1;
end

P_2_kPa
P_2 = 1000*P_2_kPa;
Z_2 = py.CoolProp.CoolProp.PropsSI('Z','P',P_2,'T',T_2,'Air');
nu_2_Po = py.CoolProp.CoolProp.PropsSI('viscosity','P',P_2,'T',T_2,'Air');

u_2 = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z_2*T_2/P_2)
Re_2 = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_2_Po*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm

%% Sanity checks

% Velocity
% D_mm = 476
% P_a = 100
% P_1 = 7000
% T_1 = 288
% T_a = 288
% Q_a_day = 7500000
% G = 0.6
% Z_1=0.95
% 
% u_1 = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z_1*T_1/P_1); % Q_b has to be converted to m3/day and D to mm

% Reynolds
% D_mm = 476
% P_a_kPa = 101
% T_a = 288
% Q_a_day = 3000000
% G = 0.6
% nu_1_po = 0.00012
% Re = 0.5134*(P_a/T_a)*(G*Q_a/(nu_1*D))

% Friction factor
% P_a_kPa = 101;
% T_a = 288;
% Q_a_day = 6000000;
% G = 0.6;
% nu_1_po = 0.00012;
% % Pipe
% L = 1000; % 1 km pipe
% D = 476; % mm
% eps = 0.03; % Absolute roughness 0.03 mm
% epsD = eps/D;
% 
% % Iterations for thet estimation of friction factor using Colebrook equation
% f = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.
% df = 10;
% count=0;
% 
% while (df > 0.0001 & count < 10) 
%     f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2);
%     df = (f_n - f)/f;
%     f = f_n;
%     count = count + 1; 
% end
% f

% Standard flow rate
% E = 0.92
% T_a = 288
% P_a_kPa = 101
% P_1_kPa = 7480
% P_2_kPa = 6000
% G = 0.6
% T_f = 293
% L_e_km = 24
% Z = 0.9
% D_mm = 288
% Q = 1.002e-2*E*(T_a/P_a_kPa)^1.02*...
%     ( (P_1_kPa^2 - P_2_kPa^2)/(G^0.961*T_f*L_e_km*Z) )^0.51*D_mm^2.53 % Panhandle B - full

% Solving flow equation to find P_2 - Panhandle B equation
% D_mm = 288
% T_a = 288
% T_f = 293
% Q_a_day = 3.5e6
% E = 0.92
% P_a_kPa = 101
% G = 0.6
% L_km = 24
% P_2 = 6000
% Z_f = 0.9
% % P_2 = sqrt( P_1^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*L_km*Z_f )
% P_1 = sqrt( P_2^2 + ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*L_km*Z_f )


%%


% Friction factors for gases
% f_W = 0.032*D^-0.333        % Waymouth
% f_PanA = 0.085*Re^-0.147    % Panhandle A - Low pressures
% f_PanB = 0.015*Re^-0.039    % Panhandle B - High pressures

% Compressiblity factor
% Z = P*v/(R*T)
% Appropriate average pressure formula P_f =2/3*( P_1+P_2-(P_1*P_2)/(P_1+P_2) )
% T_f = (T_1+T_2)/2

% Flow equations
%E = 0.95;   % Pipeline efficiency
% Waymouth
% Panhandle A
% Q = 1.002e-2*E*(T_a/P_a)^1.02*...
%     ( (P_1^2 - e^s*P_2^2)/(G^0.961*T_f*L_e*Z) )^0.51*D^2.53 % Panhandle B - full
% Q = 1.002e-2*E*(T_a/P_a)^1.02*...
%     ( (P_1^2 - P_2^2)/(G^0.961*T_f*L*Z) )^0.51*D^2.53 % Panhandle B - No Dz

% General flow equation

Q = 1.1494e-3*(T_a/P_a_kPa)*( (P_1_kPa^2-P_2_kPa^2)/(G*T_f*L*Z_f*f) )*D_mm^2.5;
P_2 = sqrt(P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*T_f*L_km*Z_f*f );
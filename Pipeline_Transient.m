%% Transient analysis - 1 km pipe with storage - Using Panhandle B flow equation

%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15+273.15; % T = 15 oC - ~280 K
P_a = 101325 ; % P = 101.325 kPa

G = 1;      % Specific gas gravity - for air G = 1

P_in = 7e6; % 7 MPa
T_1 = 5 + 273.15;
T_f = T_1;

% Flow rate
Q_a = 14*1000000/(24*3600); 
    % Conversion from mcmd (millions of cubic meters
    % per day to cubic meters per second) - real demand estimation 70 mscm/day
    % St Fergus, 14 is the proportional relative to area of one of the 3, 900
    % mm diameter pipes

% Pipe
D = 0.900; % 900 mm
dL = 1000;
L = 1000; % 1 km pipe
eps = 0.04e-3; % Absolute roughness 0.04 mm

E = 0.8;   % Pipeline efficiency
% Usually varies from 0.6 to 0.92 as a function of liquid presence and age
% "Handbook of natural gas transmission and processing"

% Tank parameters
V_tank = (pi*D^2/4)*70000; % Pipe volume

dt = 1;
%----------------------------------------------------------------------%

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;

A = pi*D^2/4;
D_mm = 1000*D;
epsD = eps/D;
f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.

P_a_kPa = P_a/1000;
L_km = L./1000;

rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');

P_1 = P_in;
P_1_kPa = P_1/1000;

P_tank(1) = P_a;
T_tank(1) = T_a;

rho_1 = py.CoolProp.CoolProp.PropsSI('D','P',P_1,'T',T_1,'Air');
nu_1 = py.CoolProp.CoolProp.PropsSI('V','P',P_1,'T',T_1,'Air');
Z_1 = py.CoolProp.CoolProp.PropsSI('Z','P',P_1,'T',T_1,'Air');
nu_1_Po = 10*nu_1;



rho_tank(1) = py.CoolProp.CoolProp.PropsSI('D','P',P_tank(1),'T',T_tank(1),'Air');
u_tank(1) = py.CoolProp.CoolProp.PropsSI('U','P',P_tank(1),'T',T_tank(1),'Air');
m_tank(1) = rho_tank(1)*V_tank;

T_2 = T_tank(1);
T_f = (T_1+T_2)/2;
P_2_kPa(1) = P_tank(1)/1000;

h_2(1) = py.CoolProp.CoolProp.PropsSI('H','P',P_tank(1),'T',T_tank(1),'Air');

P_f(1) = 2/3*( P_1_kPa + P_2_kPa(1)-(P_1_kPa*P_2_kPa(1))/(P_1_kPa+P_2_kPa(1)) );
Z_f(1) = py.CoolProp.CoolProp.PropsSI('Z','P',P_f*1000,'T',T_f,'Air');


Q(1) = 1.002e-2*E*(T_a/P_a_kPa)^1.02*...
    ( (P_1_kPa^2 - P_2_kPa(1)^2)/(G^0.961*T_f*L_km*Z_f(1)) )^0.51*D_mm^2.53; % Panhandle B - full

Q1_dot(1) = Q(1)/(24*3600); % Conversion of Q from m3/day to m3/s
m_dot(1) = Q1_dot(1)*rho_a;

m_tank(2) = m_tank(1) + m_dot(1)*dt;
u_tank(2) = ( m_tank(1)*u_tank(1) + m_dot(1)*h_2(1) ) / m_tank(2);

rho_tank(2) = m_tank(2)/V_tank;
P_tank(2) = py.CoolProp.CoolProp.PropsSI('P','T',T_tank(1),'D',rho_tank(2),'Air');



%%

%Q_a_day = Q_a*24*3600;

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



% rho(1) = py.CoolProp.CoolProp.PropsSI('D','P',P(1),'T',T_f,'Air');
% nu(1) = py.CoolProp.CoolProp.PropsSI('V','P',P(1),'T',T_f,'Air');
% Z(1) = py.CoolProp.CoolProp.PropsSI('Z','P',P(1),'T',T_f,'Air');


%%
for i=1:length(L)-1
    
    rho(i) = py.CoolProp.CoolProp.PropsSI('D','P',P(i),'T',T_f,'Air');
    nu(i) = py.CoolProp.CoolProp.PropsSI('V','P',P(i),'T',T_f,'Air');
    Z(i) = py.CoolProp.CoolProp.PropsSI('Z','P',P(i),'T',T_f,'Air');

    nu_Po(i) = 10*nu(i);

    % Velocity m/s
    u(i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(i)*T_f/P(i)); % Q_b has to be converted to m3/day and D to mm

    Re(i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm

    % Friction factor
    % Iterations for the estimation of friction factor using Colebrook equation
    % f_old = f_guess;
    % df = 10;
    % count=0;
    % while (df > 0.0001 & count < 10) 
    %     % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
    %     f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
    %     df = (f_new - f_old)/f_old;
    %     f_old = f_new;
    %     count = count + 1; 
    % end
    % f(i) = f_old;

    P_1_kPa = P(i)/1000;
    P_2_kPa = 0.98*P_1_kPa; % Initial guess
    dP = 10;
    count = 0;
    
    while (dP > 0.0001 & count < 10)
        P_f = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
        Z_f = py.CoolProp.CoolProp.PropsSI('Z','P',P_f*1000,'T',T_f,'Air');
    
        P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*T_f*(dL/1000)*Z_f );
        dP = (P(i+1) - P_2_kPa)/P_2_kPa;
        P_2_kPa = P(i+1);
        count = count + 1;
    end
    P(i+1) = 1000*P(i+1); % conversion back to Pa
    %P(i+1) = 1000*P_2_kPa;
end

figure('Color',[1 1 1])
plot(L,P)
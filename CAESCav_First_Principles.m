clear
clc

CP = py.importlib.import_module('CoolProp.CoolProp');

%--------------------------- Setup ---------------------------------




%----------------------- PROBLEM PARAMETERS -------------------------%

% Cavern dimensions (assumed cylindrical)
H = 35; % Height
D = 40; % Diameter

% Ambient conditions
P_a = 101325;
T_a = 273.15 + 25;
rho_a = CP.PropsSI('D','P',P_a,'T',T_a,'Air');

% A - Left side - Inlet 
P_in = 7e6;
T_in = 273.15 + 60;
Q_st_in = 3e5; % standard cubic meters per hour
% Q_a = Q_st_in/3600;

% m_in = rho_a*Q_a;
% m_in = 100;
m_in = 108; % Huntorf

% Initial conditions
% P_o = 101325;
% P_o = 5e6;
P_o = 4.3e6; % Huntorf
T_o = 273.15 + 25;
v_o = 0;
% v_o = v_in;

g = 9.81;
theta = 0;

%--------------------- SIMULATION PARAMETERS ------------------------%
dt = 1;
Dt = 4*3600;

A_h = pi*D^2/4;
Vol = A_h*H;

rho_in = CP.PropsSI('D','P',P_in,'T',T_in,'Air');
cp_in = CP.PropsSI('C','P',P_in,'T',T_in,'Air');

v_in = m_in/(rho_in*A_h);

%---------------------- ARRAYS INITIALIZATION ----------------------%
t = 0:dt:Dt;

n_t = Dt/dt+1;  % n of time steps

P = zeros(n_t,1);
T = zeros(n_t,1);
rho = zeros(n_t,1);
cp = zeros(n_t,1);

% Sanity check
m = zeros(n_t,1); % Total mass in pipeline
E = zeros(n_t,1); % Total energy in pipeline

%--------------------- Properties at t=0 --------------------------%
P(1) = P_o;
T(1) = T_o;
rho(1) = CP.PropsSI('D','P',P(1),'T',T(1),'Air');
cp(1) = CP.PropsSI('C','P',P(1),'T',T(1),'Air');

m(1) = rho(1)*Vol;
E(1) = rho(1)*Vol*cp(1)*T(1);

for i = 2:length(t)
    m(i) = m(i-1) + m_in*dt;
    rho(i) = m(i)/Vol;
    E(i) = E(i-1) + m_in*cp_in*T_in*dt;
    T(i) = E(i)/(m(i)*cp(i-1));
    P(i) = CP.PropsSI('P','D',rho(i),'T',T(i),'Air');
    cp(i) = CP.PropsSI('C','D',rho(i),'T',T(i),'Air');
end

figure('Color',[1 1 1])
plotyy(t,m,t,P)
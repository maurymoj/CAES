clear
clc

CP = py.importlib.import_module('CoolProp.CoolProp');

%--------------------------- Setup ---------------------------------

% CAES process
% Options:
% 'Charging';
% 'Cycle';
% 'Idle'
Process = 'Cycle';

% heat_transfer_model 
% 'Adiabatic'
% Not fully implemented yet:
% 'Steady_state'
% 'Isothermal'
heat_transfer_model = 'Adiabatic';


%----------------------- PROBLEM PARAMETERS -------------------------%
Dt = 12*3600;
Dt_charg = 7*3600;

% Cavern dimensions (assumed cylindrical)
H = 35; % Height
D = 40; % Diameter

% Ambient conditions
P_a = 101325;
T_a = 273.15 + 25;
rho_a = CP.PropsSI('D','P',P_a,'T',T_a,'Air');

% A - Left side - Inlet 
P_in = 7e6;
P_max = 7e6;
DoD = 3e6; % Depth of discharge (in terms of pressure in Pa)
P_min = P_max - DoD;

T_in = 273.15 + 60;
Q_st_in = 3e5; % standard cubic meters per hour
% Q_a = Q_st_in/3600;

% m_in = 108; % Huntorf
% m_in = rho_a*Q_a;
m_in = 100;
m_out = 100; % absolute value, sign is added later

% Ambient conditions
P_amb = 101325;
T_amb = 273.15 + 25;

% Initial conditions
% P_o = 4.3e6; % Huntorf
P_o = P_min;
T_o = T_amb;
v_o = 0;

g = 9.81;
theta = 0;

%--------------------- SIMULATION PARAMETERS ------------------------%
dt = 1;

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

if strcmp(Process,'Charging') || strcmp(Process,'Cycle')
    m_dot = m_in;
    stage = 'Charging';
elseif strcmp(Process,'Discharging')
    P_o = P_max;
    m_dot = -m_out;
    stage = 'Discharging';
end

%--------------------- Properties at t=0 --------------------------%
P(1) = P_o;
T(1) = T_o;
rho(1) = CP.PropsSI('D','P',P(1),'T',T(1),'Air');
cp(1) = CP.PropsSI('C','P',P(1),'T',T(1),'Air');

m(1) = rho(1)*Vol;
E(1) = rho(1)*Vol*cp(1)*T(1);

stage_hist = {stage};

for i = 2:length(t)
    m(i) = m(i-1) + m_dot*dt;
    rho(i) = m(i)/Vol;
    E(i) = E(i-1) + m_dot*cp_in*T_in*dt;
    T(i) = E(i)/(m(i)*cp(i-1));
    P(i) = CP.PropsSI('P','D',rho(i),'T',T(i),'Air');
    cp(i) = CP.PropsSI('C','D',rho(i),'T',T(i),'Air');
    
    if P(i) >= P_max
        m_dot = 0;
    end
    cp_in = CP.PropsSI('C','P',P(i),'T',T_in,'Air');
    
    if strcmp(Process,'Charging') 
        if strcmp(stage,'Charging') & P(i) >= P_max
            m_dot = 0;
            i_charg_end = i;
            stage = 'Idle';
            disp('Charging completed')
        end
    elseif strcmp(Process,'Cycle')
        if strcmp(stage,'Charging') & P(i) >= P_max
            m_dot = 0;
            i_charg_end = i;
            stage = 'Idle_charg';
            disp('Charging completed')
        elseif strcmp(stage,'Idle_charg') & t(i) >= Dt_charg
            m_dot = -m_out;
            stage = 'Discharging';
            disp('Discharging started')
        elseif strcmp(stage,'Discharging') & P(i) <= P_min
            m_dot = 0;
            i_disch_end = i;
            stage = 'idle_disch';
            disp('Discharging completed')
        end
        stage_hist = {stage_hist;stage};
    else
        error('Process not identified/implemented')
    end

end

figure('Color',[1 1 1])
yyaxis left
plot(t./3600,m)
ylabel('Mass [kg]')
yyaxis right
plot(t./3600,P./1e6)
ylabel('P [MPa]')
xlabel('t [h]')
grid on
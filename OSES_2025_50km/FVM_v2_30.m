clear 
% close all
% clc
disp(strcat('start time = ',string(datetime)))
% profile on
tic
CP = py.importlib.import_module('CoolProp.CoolProp'); % Simplifies coolprop calls
% Program requires Coolprop package: http://coolprop.org/coolprop/wrappers/MATLAB/index.html#matlab

%--------------------------- Setup ---------------------------------

% Type of simulation 
% 'CAESPipe'- Pipeline storage
% 'CAESCav' - Cavern
simType = 'CAESPipe';

% CAES process
% Options:
% 'Charging_L';
% 'Discharging_L';
% 'Cycle_L';
% 'Idle'
% Not implemented yet:
% 'NCycles_L';
% 'Charging_R';
% 'Discharging_R';
% 'Cycle_R';
% 'NCycles_R';
% 'Kiuchi'
Process = 'Charging_L';

% heat_transfer_model 
% 'Adiabatic'
% 'Isothermal'
% Not fully implemented yet:
% 'Kiuchi'
% 'Steady_state'
heat_transfer_model = 'Steady_state';

if strcmp(simType,'CAESPipe')
    n_nodes = 40;
    max_iter = 40;
    
    % Setting tolerance
    if strcmp(Process,'Discharging_R') || strcmp(Process,'Discharging_L')
        tol = 1e-4; % CAESPipe discharging
        tol_v = 1e-4;
    elseif strcmp(Process,'Charging_R') || strcmp(Process,'Charging_L')
        tol = 1e-4; % CAESPipe Charging
        tol_v = 1e-4;
    elseif strcmp(Process,'Cycle_L')
        tol = 1e-4;
        tol_v = 1e-4;
    elseif strcmp(Process,'Idle')
        tol = 1e-5;
        tol_v = 1e-5;
    elseif strcmp(Process,'Kiuchi')
        tol = 1e-5;
        tol_v = 1e-4;

        n_nodes = 50;
        Delta_t = 0.01*60;
    else
        error('Process not found.')
    end

elseif strcmp(simType,'CAESCav')
    n_nodes = 4;
    max_iter = 20;
    tol = 1e-8;
else
    warning('Cant identify simulation type.')
end


% Plot figures ? [0 -> no / 1 -> yes]
Save_data = 0;
Figures = 0; 
P_Corr_fig = 0;
adapt_underrelax = 0;

save_errors = 0;

% t_ramp = 0.1; % Ramp up time for inlet velocity in [s]

%----------------------- PROBLEM PARAMETERS -------------------------%

% Total simulation time
% Dt = 360; % Charging time for 5 km pipeline
% Dt = 720; % Charging + discharging time (5km 0.5m pipeline)
% Dt = 1200; % Charging L=10 km, d=0.9 m pipeline (360s = 3.44 MWh,1054 elapsed time)


if strcmp(Process,'Discharging_L') || strcmp(Process,'Discharging_R')
    % Dt = 8*3600;
    Dt = 4*3600;
elseif strcmp(Process,'Charging_L') || strcmp(Process,'Charging_R')
    % Dt = 8*3600;
    Dt = 4*3600;
elseif strcmp(Process,'Cycle_L') || strcmp(Process,'Cycle_R')
    Dt_charg = 4*3600; % Charging + idle phase duration
    Dt = 8*3600;
elseif strcmp(Process,'Kiuchi')
    Dt = 3600; % 60 min duration
else
    error('Process not identified during set up of duration.')
end

% System operational limits
P_max = 7e6;
DoD = 3e6; % Depth of discharge in terms of pressure
P_min = P_max - DoD;

% Charging rate
m_in = 100;
% T_in = 273.15+5;
T_in = 273.15 + 60;

% Discharging rate
% m_A = 417; % Huntorf
m_out = 100;

if strcmp(simType,'CAESPipe')
    % L = 213333; % Length for eq. vol with Huntorf
    % L = 100000; % Reference case
    L = 50000; % 50 km
    D = 0.9;
elseif strcmp(simType,'CAESCav')
    % L = 300; % Volume shape similar to Huntorf
    % D = 24;
    L = 50.625; % Roughly same volume as 100 km, 0.9 D pipeline
    D = 40;
end

if strcmp(Process,'Kiuchi')
    L = 5000;
    D = 0.5;
    P_max = 5e6;
    P_min = 5e6;
    P_L = 5e6;
    T_in = 273.15+25;
    warning('Kiuchi: left boundary; right boundary; boundary matrix; and boundary update under development.')
end

% Assumptions
eps = 0.04e-3; % Absolute roughness 0.04 mm

% Ambient conditions
P_amb = 101325;
T_amb = 273.15 + 25;

% !!!!!!!  GROUND TEMPERATURE SET TO AMBIENT TEMPERATURE TEST !!!!!!!!!!!
% T_ground = 273.15 + 5; % Kostowski, 2022
T_ground = T_amb;

if strcmp(Process,'Charging_L')
    % Initial conditions
    P_0 = P_min;
    T_0 = T_ground;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Inlet'; % OR M_CONST ??
    R_bound = 'Wall';
elseif strcmp(Process,'Discharging_L')
    % Initial conditions
    P_0 = P_max;
    T_0 = T_ground;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'M_const';
    R_bound = 'Wall';
elseif strcmp(Process,'Charging_R')
    % Initial conditions
    P_0 = P_min;
    T_0 = T_ground;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Wall';
    R_bound = 'Inlet';
    
elseif strcmp(Process,'Discharging_R')
    % Initial conditions
    P_0 = P_max;
    T_0 = T_ground;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Wall';
    R_bound = 'M_const';

elseif strcmp(Process,'Cycle_L')
    % Initial conditions
    P_0 = P_min;    
    T_0 = T_ground;
    v_0 = 0;
    % v_0 = v_in;

    P_0 = 4.0624189468e+06; % Final pressure after the first cycle

    L_bound = 'Inlet'; % OR M_CONST ??
    R_bound = 'Wall';

    stage = 'Charging';
    stage_hist = {stage};
elseif strcmp(Process,'Cycle_R')
elseif strcmp(Process,'NCycles_L')
elseif strcmp(Process,'NCycles_R')
elseif strcmp(Process,'Idle')
    P_0 = P_min;
    T_0 = T_amb;
    v_0 = 0;

    L_bound = 'Wall';
    R_bound = 'Wall';
elseif strcmp(Process,'Kiuchi')
    P_0 = P_min; % 5 MPa
    T_0 = 273.15 + 25; % 25 oC
    v_0 = 0;

    T_ground = T_0;

    % For Kiuchi the standard state is at 0oC and 1 atm
    rho_st = CP.PropsSI('D','T',273.15,'P',101325,'Methane');
    vol_st = 3e5/3600; % convert 300,000 sm3/h to sm3/s
    m_dot_st = rho_st*vol_st;
    L_bound = 'Wall';
    R_bound = 'Wall';

    stage = 'Idle_Kiuchi';
    stage_hist = {stage};
else
    error('Process not identified !')
end

% A - Left side boundary condition
% L_bound = 'Wall';
% L_bound = 'Inlet';
% L_bound = 'M_const';
% Not implemented yet:
% L_bound = 'P_const';
% L_bound = 'Outlet';
if strcmp(L_bound,'Inlet')
    % Q_st_in = 3e5; % standard cubic meters per hour
    % Q_a = Q_st_in/3600;
    % m_in = rho_a*Q_a;
    % m_in = 100;
    % m_L = 108; % Huntorf
    m_L = m_in;
    % v_A = ?
    P_L = P_0;
    T_L = T_in;

elseif strcmp(L_bound,'Wall')
    m_L = 0;
    v_L = 0;
elseif strcmp(L_bound,'Outlet')
    % m_A = m_R;
elseif strcmp(L_bound,'P_const')
    P_L = P_min; % 4 MPa
elseif strcmp(L_bound,'M_const')
    m_L = -m_out; % Sign indicates flow direction
else 
    error('Left boundary type not identified')
end   

% B - Right side boundary condition
% R_bound = 'Wall';
% Not implemented yet:
% R_bound = 'Outlet';
% R_bound = 'Inlet';
% R_bound = 'P_const';
% R_bound = 'M_const';
if strcmp(R_bound,'Outlet')
    % m_R = m_out;
elseif(strcmp(R_bound,'Wall'))
    m_R = 0;
    v_R = 0;
elseif(strcmp(R_bound,'Inlet'))
    m_R = -m_in; % Remember the sign to indicate the direction of flow (+ right , - left)
    P_R = P_0; % Sliding pressure at pipeline inlet
    T_R = T_ground;
    % v_B = m_B/(rho_B*A_h);
elseif strcmp(R_bound,'P_const')
    P_R = 4e6; 
    % P_B = 4.13e6; % 4.13 MPa = Huntorf https://www.sciencedirect.com/science/article/pii/S0196890420302004#s0010
elseif strcmp(R_bound,'M_const')   
    m_R = -m_out;
else 
    error('Right boundary type not identified')
end                     

k_pipe = 45.3; % [W/m K] - Wen et al. (2023) "Heat Transfer Model of Natural Gas Pipeline Based on Data Feature Extraction and First Principle Models"
k_int_coating = 0.52;
k_ext_coating = 0.4;
cp_pipe = 500; % [J/kg K]

thickness_pipe = 15.9e-3; % Still based on Wen, but in accordance with the regulations in the UK
thickness_int_coating = 0.5e-3;
thickness_ext_coating = 3e-3;

if strcmp(Process,'Kiuchi')
    fluid = 'Methane'; % Kiuchi
else
    fluid = 'Air';
end

R = 287;
g = 9.81;
theta = 0;

disp(strcat('Process = ',Process,'; Pressure = ',num2str(P_min./1e6),'-',num2str(P_max./1e6),' MPa.'))
%--------------------- SIMULATION PARAMETERS ------------------------%

if strcmp(simType,'CAESPipe')
    % dx = L/(n_nodes-1);
    dx = L/n_nodes;
    dt_max = dx/400; % 400 is representative of the sound speed, 
                     % it is higher than the maximum sound speed reached in the pipeline 
                     % to achieve a conservative value
    dt = 0.5*(Dt/ceil(Dt/dt_max)); % division of Dt in an integer number of intervals
                            % with dt smaller than dt_max

    if strcmp(Process,'Kiuchi')
        dt = Delta_t;
    end

    if dx/400 < dt % 400 upper limit for the speed of sound
        warning('dt > time needed for pressure wave to cross a node')
    end

elseif strcmp(simType,'CAESCav')
    dx = L/(n_nodes-1);
    dt_max = dx/350;
    dt = (Dt/ceil(Dt/dt_max));
else
    warning('Cant identify simulation type.')
end

Courant_number = 340*dt/dx;
if Courant_number > 1
    warning('The Courant number is greater than 1, simulation can become unstable.')
end
% disp(strcat('time step:',num2str(dt),' s'))



% Friction calculation equation
% 'Nikuradse'
% 'Colebrook'; % No significant difference and adds calculation time
friction_model = 'Nikuradse';
if strcmp(Process,'Kiuchi')
    friction_model = 'Kiuchi';
end

% Tuning
min_iter = 6; % Minimum number of iterations

% Under-relaxation (1 means no under-relaxation)
alpha = 0.5;
alpha_P = alpha;  % Pressure under-relaxation factor
alpha_v = alpha;  % Velocity under-relaxation factor
alpha_rho = alpha ;  % Density under-relaxation factor

adapt_underrelax_hist = [alpha_P alpha_v alpha_rho];

alpha_max = 0.9;
alpha_min = 0.1;

alpha_hist = [alpha_P alpha_v alpha_rho];

% Minimum absolute value of velocity to compute the relative error
v_threshold = 1e-6;

%---------------------- ARRAYS INITIALIZATION ----------------------%
t = 0:dt:Dt;

x = 0:dx:L;
x_f = [0:dx:L]';
x_n = [dx/2:dx:L]';

n_n = L/dx;       % number of nodes
% n_n = n_nodes;
n_f = n_n + 1;      % number of faces
n_t = Dt/dt+1;  % n of time steps

% Face initialization
v = zeros(n_f,n_t);
P_f = zeros(n_f,n_t);
T_f = zeros(n_f,n_t);
rho_f = zeros(n_f,n_t);
cp_f = zeros(n_f,n_t);
h_f = zeros(n_f,n_t);
s_f = zeros(n_f,n_t);
u_f = zeros(n_f,n_t);

u_sonic_f = zeros(n_f,n_t);
% drho_dP_f = zeros(n_f,n_t);
nu = zeros(n_f,n_t);


% Node initialization
v_n = zeros(n_n,n_t);
P = zeros(n_n,n_t);
T = zeros(n_n,n_t);
rho = zeros(n_n,n_t);
h = zeros(n_n,n_t);
s = zeros(n_n,n_t);
u = zeros(n_n,n_t);
cp = zeros(n_n,n_t);
u_sonic_n = zeros(n_n,n_t);
% drho_dP_n = zeros(n_n,n_t);

% equation and compare to using only mass and momentum eqs
Q = zeros(1,n_t);
Q_n = zeros(n_n,n_t);

% Sanity check
m = zeros(1,n_t); % Total mass in pipeline
E = zeros(1,n_t); % Total energy in pipeline

m_n = zeros(n_n,n_t); % Mass at each node
E_n = zeros(n_n,n_t); % Energy at each node

% Convergence speed monitoring
conv_speed = zeros(1,n_t);

%------------------ SETTING INITIAL CONDITIONS ---------------------%
% Basic calculations
epsD = eps/D;   % Relative roughness
A_h = pi*D^2/4; % Cross-section area

% Dead state parameters
P_o = P_amb;
T_o = T_amb;
h_o = CP.PropsSI('H','P',P_o,'T',T_o,fluid);
u_o = CP.PropsSI('U','P',P_o,'T',T_o,fluid);
s_o = CP.PropsSI('S','P',P_o,'T',T_o,fluid);

% Initial conditions at nodes (i -> x, j -> t)
u_0 = CP.PropsSI('U','P',P_0,'T',T_0,fluid);
s_0 = CP.PropsSI('S','P',P_0,'T',T_0,fluid);

% Thermodynamic properties over all solution space
P(:,1) = P_0*ones(n_n,1);
T(:,1) = T_0*ones(n_n,1);
rho(:,1) = CP.PropsSI('D','P',P(1,1),'T',T(1,1),fluid);
cp(:,1) = CP.PropsSI('C','P',P(1,1),'T',T(1,1),fluid); 
h(:,1) = CP.PropsSI('H','P',P(1,1),'T',T(1,1),fluid); 
s(:,1) = CP.PropsSI('S','P',P(1,1),'T',T(1,1),fluid); 
u(:,1) = CP.PropsSI('U','P',P(1,1),'T',T(1,1),fluid);
u_sonic_n(:,1) = CP.PropsSI('speed_of_sound','P',P(1,1),'T',T(1,1),fluid);

% Initial conditions at faces (i -> x, j -> t)
v(:,1) = v_0; % V profile at t=0

% UPWIND SCHEME
% velocity in faces to the nodes
v_n(:,1) = (v(1:end-1,1) >= 0).*v((1:end-1),1) ...
    +      (v(1:end-1,1) <  0).*v((2:end),1);

% Properties in nodes to faces
P_f(2:end-1,1) = (v(2:end-1,1) >= 0).*P(1:end-1,1) ...
    +            (v(2:end-1,1) <  0).*P(2:end,1);
T_f(2:end-1,1) = (v(2:end-1,1) >= 0).*T(1:end-1,1) ...
    +            (v(2:end-1,1) <  0).*T(2:end,1);
rho_f(2:end-1,1) = (v(2:end-1,1) >= 0).*rho(1:end-1,1) ...
    +            (v(2:end-1,1) <  0).*rho(2:end,1);

% P_f(end,1) = P(end,1); % Assuming v >= 0 for t=0 at x=L
% T_f(end,1) = T(end,1); 
% rho_f(end,1) = rho(end,1);


% BOUNDARY CONDITIONS 
% (Inlet, Wall, Constant pressure, and Constant flow rate)

% L boundary conditions
if strcmp(L_bound,'Inlet')
    % At "outside" node
    rho_L = CP.PropsSI('D','P',P_L,'T',T_L,fluid);
    cp_L = CP.PropsSI('C','P',P_L,'T',T_L,fluid);
    h_L = CP.PropsSI('H','P',P_L,'T',T_L,fluid);
    s_L = CP.PropsSI('S','P',P_L,'T',T_L,fluid);

    % P(1,:) = P_L;
    % T(1,:) = T_A;
    % rho(1,:) = rho_A;

    P_f(1,1) = P_L;
    T_f(1,:) = T_L;
    rho_f(1,1) = rho_L;
    
    % v_L = m_L/(rho_f(1,1)*A_h);
    % v(1,1) = m_L/(rho_f(1,1)*A_h);

elseif strcmp(L_bound,'Wall')
    v(1,:) = 0;

    P_f(1,1) = P(1,1); % Zero gradient assumption
    T_f(1,1) = T(1,1);
    rho_f(1,1) = rho(1,1);

elseif strcmp(L_bound,'P_const') 
    % Assumptions:
    % - Outflow
    % - Zero gradient for all properties
    % - Constant cross-sectional area
    P(1,:) = P_L;
    T(1,1) = T(2,1);
    rho(1,1) = rho(2,1);

    v(1,1) = v(2,1);

    P_f(1,1) = P(1,1);
    T_f(1,1) = T(1,1);
    rho_f(1,1) = rho(1,1);
elseif strcmp(L_bound,'M_const')    
    % Assumptions:
    % - Zero gradient for all properties except u-velocity
    % - Constant cross-sectional area
    % P(1,1) = P(2,1);
    % T(1,1) = T(2,1);
    % rho(1,1) = rho(2,1);

    % v_n(1,1) = m_L/(rho(1,1)*A_h);
    % 
    % % Upwind scheme
    % v(1,1) = v_n(1,1);
    % 
    % P_f(1,1) = P(1,1);
    % T_f(1,1) = T(1,1);
    % rho_f(1,1) = rho(1,1);
    
% TEST LINEAR APPROXIMATION FOR P, rho, AND T
    P(1,1) = 2*P(2,1) - P(3,1);
    T(1,1) = 2*T(2,1) - T(3,1);
    rho(1,1) = 2*rho(2,1) - rho(3,1);

    % P_f(1,1) = 2*P_f(2,1) - P_f(3,1);
    % T_f(1,1) = 2*T_f(2,1) - T_f(3,1);
    % rho_f(1,1) = 2*rho_f(2,1) - rho_f(3,1);

    P_f(1,1) = P(1,1);
    T_f(1,1) = T(1,1);
    rho_f(1,1) = rho(1,1);

    % v(1,1) = m_L/(rho_f(1,1)*A_h);
    v(1,1) = 0;
    
    % v_n(1,1) = m_L/(rho(1,1)*A_h);
    % v_n(1,1) = (v(1,1) >= 0).*v(1,1) ...
        % +      (v(1,1) <  0).*v(2,1);
% elseif strcmp(L_bound,'Kiuchi')
%     % At "outside" node
%     rho_L = CP.PropsSI('D','P',P_L,'T',T_L,fluid);
%     cp_L = CP.PropsSI('C','P',P_L,'T',T_L,fluid);
%     h_L = CP.PropsSI('H','P',P_L,'T',T_L,fluid);
%     s_L = CP.PropsSI('S','P',P_L,'T',T_L,fluid);
% 
%     P_f(1,1) = P_L;
%     T_f(1,:) = T_L;
%     rho_f(1,1) = rho_L;
% 
%     % v_L = m_L/(rho_f(1,1)*A_h);
%     % v(1,1) = m_L/(rho_f(1,1)*A_h);
elseif strcmp(L_bound,'Outlet')
    % actual implementation happens later
else
    error('Left boundary type not identified.')
end


% R boundary conditions
if strcmp(R_bound,'Inlet')
    rho_R = CP.PropsSI('D','P',P_R,'T',T_R,fluid);
    cp_R = CP.PropsSI('C','P',P_R,'T',T_R,fluid);
    h_R = CP.PropsSI('H','P',P_R,'T',T_R,fluid);
    s_R = CP.PropsSI('S','P',P_R,'T',T_R,fluid);

    % P(end,:) = P_R;
    % T(end,:) = T_R;
    % rho(end,:) = rho_R;

    P_f(end,:) = P_R;
    T_f(end,:) = T_R;
    rho_f(end,:) = rho_R;
    % At face
    % v(end,1) = m_R/(rho_f(end,1)*A_h);

elseif(strcmp(R_bound,'Wall'))
    v(end,:) = 0;

    P_f(end,1) = P(end,1); % Zero gradient assumption
    T_f(end,1) = T(end,1);
    rho_f(end,1) = rho(end,1);


elseif(strcmp(R_bound,'P_const'))    
    % Assumptions:
    % - Outflow
    % - Zero gradient for all properties
    % - Constant cross-sectional area
    P(end,:) = P_R;
    T(end,1) = T(end-1,1);
    rho(end,1) = rho(end-1,1);

    v(end,1) = v(end-1,1);

    P_f(end,1) = P(end,1);
    T_f(end,1) = T(end,1);
    rho_f(end,1) = rho(end,1);
elseif strcmp(R_bound,'M_const')    
    % Assumptions:
    % - Zero gradient for all properties except u-velocity
    % - Constant cross-sectional area
    P(end,1) = P(end-1,1);
    T(end,1) = T(end-1,1);
    rho(end,1) = rho(end-1,1);

    v_n(end,1) = m_R/(rho(end,1)*A_h);

    % Upwind scheme
    v(end,1) = v_n(end,1);
    
    P_f(end,1) = P(end,1);
    T_f(end,1) = T(end,1);
    rho_f(end,1) = rho(end,1);    
elseif strcmp(R_bound,'Kiuchi')    
    % Wall
    v(end,:) = 0;

    P_f(end,1) = P(end,1); % Zero gradient assumption
    T_f(end,1) = T(end,1);
    rho_f(end,1) = rho(end,1);

elseif strcmp(R_bound,'Outlet')
    % actual implementation happens later
else
    error('Right boundary type not identified.')
end



% Boundary conditions (Outlet - depend on the other side of the pipeline)

% L outlet boundary conditions
if strcmp(L_bound,'Outlet') 
    % Assumptions:
    % - Zero gradient for all properties except u-velocity
    % - Velocity at neighbouting node is corrected to ensure conservation of
    % mass with the inlet
    % - Constant cross-sectional area
    P(1,1) = P(2,1);
    T(1,1) = T(2,1);
    rho(1,1) = rho(2,1);

    v_n(1,1) = (rho(end,1)*v_n(end,1))/rho(1,1); % Velocity correction

    % Upwind scheme
    v(1,1) = v_n(1,1);
    
    P_f(1,1) = P(1,1);
    T_f(1,1) = T(1,1);
    rho_f(1,1) = rho(1,1);
    
    % m_A = m_R;
end  

% R outlet boundary conditions
if strcmp(R_bound,'Outlet')
    % Assumptions:
    % - Zero gradient for all properties except u-velocity
    % - Velocity at neighbouting node is corrected to ensure conservation of
    % mass with the inlet
    % - Constant cross-sectional area
    P(end,1) = P(end-1,1);
    T(end,1) = T(end-1,1);
    rho(end,1) = rho(end-1,1);

    v_n(end,1) = (rho(1,1)*v_n(1,1))/rho(end,1); % Velocity correction

    % Upwind scheme
    v(end,1) = v_n(end,1);
    
    P_f(end,1) = P(end,1);
    T_f(end,1) = T(end,1);
    rho_f(end,1) = rho(end,1);
end 

% Thermodynamic properties
cp_f(1,1) = CP.PropsSI('C','P',P_f(1,1),'D',rho_f(1,1),fluid);
h_f(1,1) = CP.PropsSI('H','P',P_f(1,1),'D',rho_f(1,1),fluid);
s_f(1,1) = CP.PropsSI('S','P',P_f(1,1),'D',rho_f(1,1),fluid);
u_sonic_f(1,1) = CP.PropsSI('speed_of_sound','P',P_f(1,1),'D',rho_f(1,1),fluid);

cp_f(end,1) = CP.PropsSI('C','P',P_f(end,1),'D',rho_f(end,1),fluid);
h_f(end,1) = CP.PropsSI('H','P',P_f(end,1),'D',rho_f(end,1),fluid);
s_f(end,1) = CP.PropsSI('S','P',P_f(end,1),'D',rho_f(end,1),fluid);
u_sonic_f(end,1) = CP.PropsSI('speed_of_sound','P',P_f(end,1),'D',rho_f(end,1),fluid);

for i = 1:n_n  % PROPERTIES FROM P AND RHO          
    % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),fluid);
    cp(i,1) = CP.PropsSI('C','P',P(i,1),'D',rho(i,1),fluid);
    h(i,1) = CP.PropsSI('H','P',P(i,1),'D',rho(i,1),fluid);
    s(i,1) = CP.PropsSI('S','P',P(i,1),'D',rho(i,1),fluid);
    u(i,1) = CP.PropsSI('U','P',P(i,1),'D',rho(i,1),fluid);
    cp_f(i,1) = CP.PropsSI('C','P',P_f(i,1),'D',rho_f(i,1),fluid);
    h_f(i,1) = CP.PropsSI('H','P',P_f(i,1),'D',rho_f(i,1),fluid);
    s_f(i,1) = CP.PropsSI('S','P',P_f(i,1),'D',rho_f(i,1),fluid);
    u_sonic_f(i,1) = CP.PropsSI('speed_of_sound','P',P_f(i,1),'D',rho_f(i,1),fluid);
end

m(1) = sum(rho(:,1)*A_h*dx);
E(1) = sum(rho(:,1)*A_h*dx.*u(:,1));

m_n(:,1) = rho(:,1)*A_h*dx;
E_n(:,1) = rho(:,1)*A_h*dx.*u(:,1);

n_iters = zeros(n_t,1);

% Friction factor based on Nikuradse
f_guess = (2*log10(1/epsD)+1.14)^(-2);

if P_Corr_fig    
    resFig = figure;
end

if save_errors || P_Corr_fig
    error_hist = [];
    error_hist_v = [];
    error_hist_rho = [];
    % error_hist2 = zeros(n_n,n_t,max_iter);
end
bound_hist = [string(L_bound) string(R_bound)];

P_is_conv = zeros(n_t,1);
v_is_conv = zeros(n_t,1);

for j=2:n_t
    % Initial guess for next time step is the same props as the previous t step
    P(:,j) = P(:,j-1);
    rho(:,j) = rho(:,j-1);
    T(:,j) = T(:,j-1);

    v(:,j) = v(:,j-1);
    % v(2:end-1,j) = v(2:end-1,j-1);

    P_f(:,j) = P_f(:,j-1);
    rho_f(:,j) = rho_f(:,j-1);
    T_f(:,j) = T_f(:,j-1);

    v_n(:,j) = v_n(:,j-1);

    P_corr = zeros(n_n,1);
    rho_corr = zeros(n_n,1);
    v_corr = zeros(n_f,1);
    
    v_star = v(:,j);
    n_iters(j) = 0;
    error_P = 10;
    error_rho = 10;
    error_v = 10*ones(n_f,1);

    prev_error_P = error_P; 
    prev_error_v = error_v; 
    prev_error_rho = error_rho; 

    % Momentum and mass balance loop
    % while count(j) < max_iter && max(abs(error_P)) > tol
    while (n_iters(j) < max_iter && (max(abs(error_P)) > tol || max(abs(error_v)) > tol_v)) ...
           || (n_iters(j) <= min_iter)
        % Under-relaxed corrections
        P(:,j) = P(:,j) + alpha_P*P_corr;
        rho(:,j) = alpha_rho*(rho(:,j) + rho_corr) ...
            + (1-alpha_rho)*rho(:,j);
        v(:,j) = alpha_v*(v_star + v_corr) ...
            + (1-alpha_v)*v_star;

        if strcmp(heat_transfer_model,'Isothermal')
            for i=1:n_n
                rho(i,j) = CP.PropsSI('D','P',P(i,j),'T',T_ground,'air');
                % 
            end
        end

        % v(1:end-1,j) = alpha_v*(v_star(1:end-1) + v_corr(1:end-1)) ...
        % + (1-alpha_v)*v_star(1:end-1);
        
        alpha_hist = [alpha_hist;alpha_P alpha_v alpha_rho];

        % Boundary conditions
        if strcmp(L_bound,'Inlet')
            % P(1,:) = P_L;
            % T(1,:) = T_L;
            % rho(1,:) = rho_L;
            
            % SLIDING PRESSURE
            P_f(1,j) = P(1,j); % Assumption of constant pressure between 
                               % face and first node
            T_f(1,j) = T_L; 
            % T_f(1,j) = T(1,j); 
            rho_f(1,j) = CP.PropsSI('D','P',P_f(1,j),'T',T_f(1,j),fluid);
            v(1,j) = m_L/(rho_f(1,j)*A_h); % Sliding pressure
            
            v_n(1,j) = m_L/(rho(1,j)*A_h);

        elseif strcmp(L_bound,'Wall')
            v(1,j) = 0;
            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);

            v_n(1,j) = (v(1,j) >= 0).*v(1,j) ...
            +      (v(1,j) <  0).*v(2,j);
        elseif strcmp(L_bound,'P_const') % Kiuchi
            
            
            % warning('Kiuchi boundary conditions under implementation.')
            
            P(1,j) = P_L;
            T(1,j) = T(2,j);
            rho(1,j) = rho(2,j);

            v(1,j) = v(2,j);
        
            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);

            v_n(1,j) = v(1,j);

        elseif strcmp(L_bound,'M_const') 
            % Variation of outlet bound. condition, 
            % with known mass flow rate
            % Assumptions:
            % - Zero gradient for all properties except u-velocity
            % - Constant cross-sectional area
            % P(1,j) = P(2,j);
            % T(1,j) = T(2,j);
            % rho(1,j) = rho(2,j);
            % 
            % v_n(1,j) = m_L/(rho(1,j)*A_h);
            % 
            % % Upwind scheme
            % v(1,j) = v_n(1,j); % zero gradient assumption
            % 
            % P_f(1,j) = P(1,j);
            % T_f(1,j) = T(1,j);
            % rho_f(1,j) = rho(1,j);

            % TEST LINEAR APPROXIMATION FOR P, rho, AND T
            % !!! ASSUMING FLOW ALWAYS TO THE LEFT !!!

            P(1,j) = 2*P(2,j) - P(3,j);
            T(1,j) = 2*T(2,j) - T(3,j);
            rho(1,j) = 2*rho(2,j) - rho(3,j);
        
            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);
             
            % P_f(1,j) = 1.5*P(1,j) - 0.5*P(2,j);
            % T_f(1,j) = 2*T(1,j) - T(2,j);
            % rho_f(1,j) = 2*rho(1,j) - rho(2,j);
            
            % UNDER DEVELOPMENT
            % v(1,j) = min(1 , t(j)/t_ramp)*m_L/(rho_f(1,j)*A_h); % Ramps
            % mass flow rate from 0 up to max value based on time t_ramp
            v(1,j) = m_L/(rho_f(1,j)*A_h);

            % Node velocity differencing scheme 
            % Zero-gradient assumption
            % v_n(1,j) = v(1,j);
            % Upwind scheme
            v_n(1,j) = (v(1,j) >= 0).*v(1,j) ...
                +      (v(1,j) <  0).*v(2,j);
            % Central scheme
            % v_n(1,j) = (v(1,j) + v(2,j))/2;
            % Hybrid scheme
            % nu_f = CP.PropsSI('viscosity','P',P_f(1,j),'D',rho_f(1,j),fluid);
            % Pe = abs(rho_f(1,j)*v(1,j)*dx)/nu_f;
            % if Pe > 2  % High convection, use more upwind
            %     beta = 0.2;  % More upwind
            % elseif Pe < 0.5  % Low convection, use more central
            %     beta = 0.8;  % More central
            % else
            %     beta = 0.5;  % Balanced
            % end
            % 
            % v_n(1,j) = (1 - beta) * (v(1,j) >= 0) * v(1,j) + ...
            %            (1 - beta) * (v(1,j) <  0) * v(2,j) + ...
            %                  beta * (v(1,j) + v(2,j)) / 2;
            % simplified QUICK scheme
            % v_n(1,j) = v(1,j) + (1/8) * (3*v(2,j) - 3*v(1,j) + v(3,j) - v(1,j));
        end

        % R boundary
        if strcmp(R_bound,'Inlet')
            % P(end,:) = P_R;
            % T(end,:) = T_R;
            % rho(end,:) = rho_R;

            P_f(end,j) = P(end,j); % Assumption of constant pressure between 
                               % face and first node
            T_f(end,j) = T_R; 
            % T_f(end,j) = T(end,j); 
            rho_f(end,j) = CP.PropsSI('D','P',P_f(end,j),'T',T_f(end,j),fluid);
            v(end,j) = m_R/(rho_f(end,j)*A_h); % Sliding pressure test
            
            v_n(end,j) = m_R/(rho(end,j)*A_h);
        elseif(strcmp(R_bound,'Wall'))
            % Set by upwind scheme
            v(end,j) = 0;

            P_f(end,j) = P(end,j); % Zero gradient assumption
            T_f(end,j) = T(end,j);
            rho_f(end,j) = rho(end,j);

            % % v_n(end,j) = v(end-1,j);
            % v_n(end,j) = (v(end-1,j) >= 0).*v(end-1,j) ...
            % +      (v(end-1,j) <  0).*v(end,j); 
            % DOUBLE-CHECK !!!!!!

        elseif strcmp(R_bound,'P_const')
            % P(end,j) = P(end-1,j);
            P(end,j) = P_R;
            T(end,j) = T(end-1,j);
            rho(end,j) = rho(end-1,j);
        
            v(end,j) = v(end-1,j);
        
            P_f(end,j) = P(end,j);
            T_f(end,j) = T(end,j);
            rho_f(end,j) = rho(end,j); 
            
            v_n(end,j) = v(end,j);
        elseif strcmp(R_bound,'M_const')  
            % Assumptions:
            % - Zero gradient for all properties except u-velocity
            % - Constant cross-sectional area
            P(end,j) = P(end-1,j);
            T(end,j) = T(end-1,j);
            rho(end,j) = rho(end-1,j);
        
            v_n(end,j) = m_R/(rho(end,j)*A_h);
        
            % Upwind scheme
            v(end,j) = v_n(end,j);
            
            P_f(end,j) = P(end,j);
            T_f(end,j) = T(end,j);
            rho_f(end,j) = rho(end,j); 
        elseif strcmp(R_bound, 'Vol_const')


            % warning('Kiuchi right boundary under implementation.')
            

            % Zero gradient assumption
            P(end,j) = P(end-1,j);
            T(end,j) = T(end-1,j);
            rho(end,j) = rho(end-1,j);
        
            P_f(end,j) = P(end,j);
            T_f(end,j) = T(end,j);
            rho_f(end,j) = rho(end,j);

            % Calculation of mass flow rate
            % m_dot_st = rho_st*vol_st; % Standard mass flow rate
            v(end,j) = m_dot_st/(rho_f(end,j)*A_h);
            
            v_n(end,j) = v(end,j);
            % v_n(end,j) = (rho(1,j)*v_n(1,j))/rho(end,j); % Velocity correction            

        end 

        if strcmp(L_bound,'Outlet')
            P(1,j) = P(2,j);
            T(1,j) = T(2,j);
            rho(1,j) = rho(2,j);

            v_n(1,j) = rho(end,j)*v(end,j)/rho(1,j);
            v(1,j) = v_n(1,j);

            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);
        elseif strcmp(R_bound,'Outlet')
            P(end,j) = P(end-1,j);
            T(end,j) = T(end-1,j);
            rho(end,j) = rho(end-1,j);
        
            v_n(end,j) = (rho(1,j)*v_n(1,j))/rho(end,j); % Velocity correction
            v(end,j) = v_n(end,j); % ENSURING CONSERVATION OF MASS

            P_f(end,j) = P(end,j);
            T_f(end,j) = T(end,j);
            rho_f(end,j) = rho(end,j);
        end

        % Upwind scheme
        % v_n(:,j) = (v(1:end-1,j) >= 0).*v((1:end-1),j) ...
        %     +      (v(1:end-1,j) <  0).*v((2:end),j);
        v_n(2:end,j) = (v(2:end-1,j) >= 0).*v((2:end-1),j) ...
            +      (v(2:end-1,j) <  0).*v((3:end),j);
        % QUICK + Upwind scheme
        % v_n(2:5,j) = v(2:5,j) + (1/8) * (3*v(3:6,j) - 3*v(2:5,j) + v(4:7,j) - v(1:4,j));
        % v_n(6:end,j) = (v(6:end-1,j) >= 0).*v(6:end-1,j) ...
        %     +      (v(6:end-1,j) <  0).*v(7:end,j);
        % QUICK scheme
        % v_n(i)       = v_f(i)       + (1/8) * (3*v_f(i+1)   - 3*v_f(i)       + v_f(i+2)     - v_f(i-1))
        % v_n(2:n_n-1,j) = v(2:n_n-1,j) + (1/8) * (3*v(3:n_n,j) - 3*v(2:n_n-1,j) + v(4:n_n+1,j) - v(1:n_n-2,j));
        

        % Properties in nodes to faces
        P_f(2:end-1,j) = (v(2:end-1,j) >= 0).*P(1:end-1,j) ...
            +            (v(2:end-1,j) <  0).*P(2:end,j);
        rho_f(2:end-1,j) = (v(2:end-1,j) >= 0).*rho(1:end-1,j) ...
            +            (v(2:end-1,j) <  0).*rho(2:end,j);


        % Properties
        % u_sonic_f = zeros(n_f,1);
        drho_dP_f = zeros(n_f,1);
        nu = zeros(n_f,1);

        % u_sonic_n = zeros(n_n,1);
        drho_dP_n = zeros(n_n,1);

        for i = 1:n_n  % PROPERTIES FROM P AND RHO
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),fluid);

            % u_sonic_f(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),fluid);
            % if strcmp(friction_model,'Colebrook')
            %     nu = CP.PropsSI('viscosity','P',P_f(i,j),'D',rho_f(i,j),fluid);
            % end
            % 
            % drho_dP_f(i) = 1/(u_sonic_f(i)^2);
            % 
            % u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),fluid);
            % drho_dP_n(i) = 1/(u_sonic_n(i)^2);

            % u_sonic_f(i,j) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),fluid);
            if strcmp(friction_model,'Colebrook')
                nu(i) = CP.PropsSI('viscosity','P',P_f(i,j),'D',rho_f(i,j),fluid);
            end

            drho_dP_f(i) = 1/(u_sonic_f(i,j-1)^2);

            % u_sonic_n(i,j) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),fluid);
            drho_dP_n(i) = 1/(u_sonic_n(i,j-1)^2);

        end
        
        % u_sonic_f(end,j) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'D',rho_f(end,j),fluid);
        drho_dP_f(end) = 1/(u_sonic_f(end,j-1)^2);
        if strcmp(friction_model,'Colebrook')
            nu(end) = CP.PropsSI('viscosity','P',P_f(end,j),'D',rho_f(end,j),fluid);
        end

        %---------- v* calculation - Momentum control volume ---------------------%
        
        % Friction factor based on Nikuradse
        if strcmp(friction_model,'Nikuradse')
            f = (2*log10(1/epsD)+1.14)^(-2)*ones(n_f,1);
        elseif strcmp(friction_model,'Colebrook')
            Re = rho_f(:,j).*abs(v(:,j))*D./nu(:);
            f = zeros(n_f,1);
    
            for i=1:n_f
                f_old = f_guess;
                df = 10;
                count_f = 0;
                while (df > 0.0001 & count_f < 20) 
                    f_new = (-2*log10(epsD/3.7 + 2.51/(Re(i)*f_old^0.5)))^(-2); % Original Colebrook-White equation
                    % f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
                    df = abs((f_new - f_old)/f_old);
                    f_old = f_new;
                    count_f = count_f + 1; 
                end
                f(i) = f_old;
            end
        elseif strcmp(friction_model,'Kiuchi')
            f = 0.008*ones(n_f,1);
        else
            error('Friction model not identified')
        end
        
        

        % Matrix/vector initialization (only need to solve for the inner faces -
        % boundary faces are dealt with the boundary conditions)
        a_M = zeros(n_f);
        b_M = zeros(n_f,1);
        d_M = zeros(n_f,1);
        B_M = zeros(n_f,1);
        

        a_M(1,1) = 1e10; % Value for v at the first momentum volume is 
                         % known from boundary conditions
        B_M(1) = 1e10*v(1,j);
        %!!!!!!!!!!!!!!!!!!! IN DEVELOPMENT !!!!!!!!!!!!!!!!!!!!!!!!!
        % if strcmp(L_bound,'P_const') 
        %     a_M(1,1) = 1; % Value for v at the first momentum volume is 
        %               % known from boundary conditions
        %     B_M(1) = v(1,j);
        % end
        %!!!!!!!!!!!!!!!!!!! IN DEVELOPMENT !!!!!!!!!!!!!!!!!!!!!!!!!

        % Right boundary
        a_M(n_f,n_f) = 1; % Value for v at the last momentum volume
        B_M(n_f) = v(n_f,j);
        
        % 
        a_M(2,3) = -max( 0, -rho(2,j)*v_n(2,j)/dx);    % a_C
        a_M(2,2) = rho_f(2,j)/dt + f(2)*rho_f(2,j)*abs(v(2,j))/(2*D) ...
            + (rho(2,j)*v_n(2,j) - rho(1,j)*v_n(1,j))/dx ...
            - (a_M(2,3) - max(rho(1,j)*v_n(1,j)/dx,  0) );% a_B

        d_M(2) = -1/dx;

        b_M(2) = (rho_f(2,j-1)*v(2,j-1)/dt)...
            - rho_f(2,j)*g*sind(theta)...
            + max(rho(1,j)*v_n(1,j)/dx,  0)*v_n(1,j); 
                                % SHOULD WE USE V(1) OR V_N(1) HERE ?

        B_M(2) = d_M(2)*(P(2,j)-P(1,j)) ...
            + b_M(2);



        % if strcmp(L_bound,'M_const') % Assumed outlet
        %     a(1,2) = rho(2,j)*v_n(2,j)/dx;    % a_C
        %     a(1,1) = rho_f(2,j)/dt + f(1)*rho_f(2,j)*abs(v(2,j))/(2*D) ...
        %         - rho(1,j)*v_n(1,j)/dx;% a_B
        % 
        %     d(1) = -1/dx;
        %     b(1) = (rho_f(2,j-1)*v(2,j-1)/dt)...
        %         - rho_f(2,j)*g*sind(theta);
        % 
        %     B(1) = d(1)*(P(2,j)-P(1,j)) ...
        %         + b(1);  
        % end

        % IMPLEMENT M_CONST BOUNDARY FOR R_bound




        a_M(n_f-1,n_f-2) = -max(rho(n_n-1,j)*v_n(n_n-1,j)/dx, 0); % a_N-2

        a_M(n_f-1,n_f-1) = ...                                  % a_N-1
            rho_f(n_f-1,j)/dt + f(n_f-1)*rho_f(n_f-1,j)*abs(v(n_f-1,j))/(2*D)...
            +(rho(n_n,j)*v_n(n_n,j) - rho(n_n-1,j)*v_n(n_n-1,j))/dx...
            - (a_M(n_f-1,n_f-2) - max(0, -rho(n_n,j)*v_n(n_n,j)/dx) );
        
        d_M(n_f-1) = -1/dx;
        b_M(n_f-1) = rho_f(n_f-1,j-1)*v(n_f-1,j-1)/dt...
            -rho_f(n_f-1,j)*g*sind(theta)...
            + max(0, -rho(n_n,j)*v_n(n_n,j)/dx)*v_n(n_n,j);
        
        B_M(n_f-1) = d_M(n_f-1)*(P(n_n,j)-P(n_n-1,j)) + b_M(n_f-1);
        

        

        for i=3:n_f-2

            a_M(i,i-1) = -max(rho(i-1,j)*v_n(i-1,j)/dx,                     0);
            a_M(i,i+1) = -max(                       0, -rho(i,j)*v_n(i,j)/dx);
            a_M(i,i)   = rho_f(i,j)/dt + f(i)*rho_f(i,j)*abs(v(i,j))/(2*D) ...
                + ( rho(i,j)*v_n(i,j) - rho(i-1,j)*v_n(i-1,j) )/dx...
                - ( a_M(i,i-1) + a_M(i,i+1) );
            d_M(i) = -1/dx;
            b_M(i) = (rho_f(i,j-1)*v(i,j-1)/dt)-rho_f(i,j)*g*sind(theta);
            
            B_M(i) = d_M(i)*(P(i,j)-P(i-1,j)) + b_M(i);
        end
        
        v_star = linsolve(a_M,B_M);
                
        % if strcmp(L_bound,'M_const')
        %     v_star(1) = v(1,j);
        % end
        % if strcmp(L_bound,'Outlet')
        %     v_star(1) = v_star(2);
        % elseif strcmp(R_bound,'Outlet')
        %     v_star(end) = v_star(end-1);
        % end


        %--------------- Pressure correction P' calculations ---------------------%
        a_C = zeros(n_n);
        B_C = zeros(n_n,1);
        
        % left node
        a_C(1,2) = (rho_f(2,j)*d_M(2)/a_M(2,2))/dx ...
            - max(0, -drho_dP_f(2)*v_star(2)/dx); % A_2
        a_C(1,1) = drho_dP_n(1)/dt + drho_dP_f(2)*v_star(2)/dx ...
               - a_C(1,2); % A_1
        
        B_C(1) = (rho(1,j-1)-rho(1,j))/dt ...
            + (rho_f(1,j)*v_star(1) - rho_f(2,j)*v_star(2))/dx;
        
        if strcmp(L_bound,'P_const')
            a_C(1,1) = 1;
            a_C(1,2) = 0;

            B_C(1) = 0;
        elseif strcmp(L_bound,'Outlet')
            %???
            a_C(1,2) = -drho_dP_f(2)*v_star(2)/dx + rho_f(2,j)*(d_M(1)/a_M(1,1))/dx; % A_2
            a_C(1,1) = drho_dP_n(1)/dt ...
                     - (rho_f(2,j)*d_M(1)/a_M(1,1))/dx ...
                     - drho_dP_f(1)*v_star(1)/dx ; % A_1
        elseif strcmp(L_bound,'Kiuchi')
            error('Kiuchi left boundary matrix coefficients not set yet.')
            % BOUNDARY CHANGES WITH TIME
            % Constant pressure
            a_C(1,1) = 1;
            a_C(1,2) = 0;

            B_C(1) = 0;
        end
        
        
        % right node
        a_C(n_n,n_n-1) = (rho_f(n_f-1,j)*d_M(n_f-1)/a_M(n_f-1,n_f-1))/dx...
                       - max(drho_dP_f(n_f-1)*v_star(n_f-1)/dx, 0);
        a_C(n_n,n_n) = drho_dP_n(n_n)/dt ...
                     - drho_dP_f(n_f-1)*v_star(n_f-1)/dx... 
                     - a_C(n_n,n_n-1); % v_star(end) = 0
        
        B_C(n_n) = (rho(n_n,j-1)-rho(n_n,j))/dt ...
            + (rho_f(n_f-1,j)*v_star(n_f-1) - rho_f(n_f,j)*v_star(n_f))/dx;
        
        

        if strcmp(R_bound,'P_const')
            a_C(n_n,n_n-1) = 0;
            a_C(n_n,n_n) = 1;

            B_C(n_n) = 0;
        elseif strcmp(R_bound,'Outlet')
            %???
            a_C(n_n,n_n-1) = (rho_f(N_f-1,j)*d_M(end)/a_M(end,end))/dx...
                           - drho_dP_f(N_f-1)*v_star(N_f-1)/dx;
            a_C(n_n,n_n) = drho_dP_n(n_n)/dt ...
                         - (rho_f(N_f-1,j)*d_M(end)/a_M(end,end))/dx ...
                         + drho_dP_f(N_f)*v_star(N_f)/dx;
        elseif strcmp(R_bound,'Kiuchi')
            error('Kiuchi right boundary matrix coefficients not set yet.')
            % Outlet
            a_C(n_n,n_n-1) = (rho_f(N_f-1,j)*d_M(end)/a_M(end,end))/dx...
                           - drho_dP_f(N_f-1)*v_star(N_f-1)/dx;
            a_C(n_n,n_n) = drho_dP_n(n_n)/dt ...
                         - (rho_f(N_f-1,j)*d_M(end)/a_M(end,end))/dx ...
                         + drho_dP_f(N_f)*v_star(N_f)/dx;
        end



        for i=2:n_n-1
            a_C(i,i-1) = (rho_f(i,j)*d_M(i)/a_M(i,i))/dx ...
                       - max(drho_dP_f(i)*v_star(i)/dx,           0);
            a_C(i,i+1) = (rho_f(i+1,j)*d_M(i+1)/a_M(i+1,i+1))/dx ...
                       - max( 0, -drho_dP_f(i+1)*v_star(i+1)/dx);
            a_C(i,i)   = drho_dP_n(i)/dt + (drho_dP_f(i+1)*v_star(i+1)/dx ...
                       - drho_dP_f(i)*v_star(i)/dx) - ( a_C(i,i-1) + a_C(i,i+1) );
            
            B_C(i) = (rho(i,j-1)-rho(i,j))/dt ...
                    + (rho_f(i,j)*v_star(i) - rho_f(i+1,j)*v_star(i+1))/dx;
        end

        P_corr = linsolve(a_C,B_C);
        
        if strcmp(L_bound,'Inlet')
            % P_corr(1) = 0;
        elseif strcmp(R_bound,'M_const')
            P_corr(end) = 0;
        elseif strcmp(L_bound,'Kiuchi')
            error('Kiuchi pressure correction at boundary not set yet.')
        end

        rho_corr = drho_dP_n.*P_corr;

        a_diag = diag(a_M);


        v_corr = zeros(n_f,1);
        % v_corr(1)         = d(1)./a_i(1).*P_corr(1);
        v_corr(2:end-1)   = d_M(2:end-1)./a_diag(2:end-1).*(P_corr(2:end)-P_corr(1:end-1));

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Inlet boundary condition
        % v_corr(1)         = 0; 
        % Outlet boundary condition
        % v_corr(end)       = 0;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        % v_corr(end) = % Velocity corrected based on mass balance (m_in =
        % m_out)
        error_P = (P_corr./P(:,j));
        error_rho = (rho_corr./rho(:,j));
        
        v_mask = abs(v(:,j)) > v_threshold; % mask to use absolute values when v is too small or 0
        error_v(v_mask) = (v_corr(v_mask)./v(v_mask,j));
        error_v(~v_mask) = abs(v_corr(~v_mask));
        
        if save_errors || P_Corr_fig
            error_hist = [error_hist error_P];
            % error_hist2(:,j,n_iters(j)+1) = error_P;
        
            error_hist_rho = [error_hist_rho error_rho];
            error_hist_v = [error_hist_v error_v];
        end

        if rem(n_iters(j),5) == 0 && P_Corr_fig
            figure(resFig)
            plot(mean(abs(error_hist)))
            % ylim([0,1e-4])
            drawnow % limitrate
        end


        if any(isnan(P_corr)) || any(isinf(P_corr))
            warning('Divergence detected: Solution contains NaN or Inf values');
            break;  % or take corrective action
        end
        
        % if count(j) > min_iter && std(P_corr(end-min_iter+1:end)) < std(P_corr(end-min_iter*2+1:end-min_iter))
        %     warning('Possible oscillatory behavior detected');
        %     % Adjust under-relaxation factors or time step
        % end

        % Adaptable under-relaxation factor based on error relative to
        % tolerance
        if adapt_underrelax
            if max(abs(error_P)) > tol * 10  % if error is much larger than tolerance
                alpha_P = max(alpha_P /2 , alpha_min);  % reduce relaxation factor, but not below a minimum
            elseif max(abs(error_P)) < prev_error_P / 2
                alpha_P = min(alpha_P *2, alpha_max);  % increase relaxation factor, but not above 1
            end
    
            if max(abs(error_v)) > tol_v * 10  % if error is much larger than tolerance
                alpha_v = max(alpha_v /2, alpha_min);  % reduce relaxation factor, but not below a minimum
            elseif max(abs(error_v)) < prev_error_v / 2
                alpha_v = min(alpha_v *2, alpha_max);  % increase relaxation factor, but not above 1
            end
    
            if max(abs(error_rho)) > tol * 10  % if error is much larger than tolerance
                alpha_rho = max(alpha_rho /2, alpha_min);  % reduce relaxation factor, but not below a minimum
            elseif max(abs(error_rho)) < prev_error_rho / 2
                alpha_rho = min(alpha_rho *2, alpha_max);  % increase relaxation factor, but not above a maximum
            end

            prev_error_P = max(abs(error_P));
            prev_error_v = max(abs(error_v));
            prev_error_rho = max(abs(error_rho));
            adapt_underrelax_hist = [adapt_underrelax_hist;alpha_P alpha_v alpha_rho];
        end
        
        prev_error_P = max(abs(error_P));
        
        n_iters(j) = n_iters(j)+1;
    end
    
    if n_iters(j) >= 0.75*max_iter
        conv_speed(j) = -1; % Convergence is slow, decrease time step
    elseif n_iters(j) <= 0.25*max_iter
        conv_speed(j) = 1; % Convergence is fast, increase time step 
    end

    if max(abs(error_P)) >= tol 
        P_is_conv(j) = 0;
    else
        P_is_conv(j) = 1;
    end

    if max(abs(error_v)) >= tol_v 
        v_is_conv(j) = 0;
    else
        v_is_conv(j) = 1;
    end

    % alpha_P  = alpha;
    % alpha_v  = alpha + 0.4;
    % alpha_rho = alpha;

    % ENERGY BALANCE
    % ASSUMING v>0
    a_T = zeros(n_n,n_n);
    b_T = zeros(n_n,1);

    if strcmp(heat_transfer_model,'Adiabatic')
        Q_n(:,j) = 0;
        Q(j) = 0;
    elseif strcmp(heat_transfer_model,'Isothermal')

    elseif strcmp(heat_transfer_model,'Sukhov')

    elseif strcmp(heat_transfer_model,'Steady_state')
        % neglecting effect of coating
        D_ext = D + thickness_pipe; 
        U = dx*k_pipe/(D*log(D_ext/D));
        Q_n(:,j) = -4*U/D*(T(:,j-1) - T_ground);
        Q(j) = sum(Q_n(:,j));
    elseif strcmp(heat_transfer_model,'Transient')

    elseif strcmp(heat_transfer_model, 'Kiuchi')
        % U = 2.84; % [W/m2K] Chaczykowski, 2010
        U = 0; % Adiabatic
        Q_n(:,j) = -4*U/D*(T(:,j-1) - T_ground);
        Q(j) = sum(Q_n(:,j));
    else
        error('Heat transfer model not identified')
    end      
    
    if strcmp(heat_transfer_model,'Isothermal')
        T(:,j) = T_ground;
    else
        a_T(1,2) = -max(0, -rho_f(2,j)*v(2,j)*cp_f(2,j-1)/dx);    
    
        a_T(1,1) = rho(1,j)*cp(1,j-1)/dt ...
            + (rho_f(2,j)*v(2,j)*cp_f(2,j-1) - rho_f(1,j)*v(1,j)*cp_f(1,j-1))/dx...
            - (a_T(1,2) - max(rho_f(1,j)*v(1,j)*cp_f(1,j-1)/dx, 0));
        
        b_T(1) = Q_n(1,j) + (P(1,j)-P(1,j-1))/dt ...
                + v_n(1,j)*(P_f(2,j) - P_f(1,j))/dx ...
                + f(1)*rho(1,j)*abs(v_n(1,j))^3/(2*D) ...
                + rho(1,j-1)*cp(1,j-1)*T(1,j-1)/dt...
                + max(rho_f(1,j)*v(1,j)*cp_f(1,j-1)/dx, 0)*T_f(1,j);
    
    
        a_T(n_n,n_n-1) = -max(rho_f(n_f-1,j)*v(n_f-1,j)*cp_f(n_f-1,j-1)/dx, 0);
    
        a_T(n_n,n_n) = rho(n_n,j)*cp(n_n,j-1)/dt ...
            + (rho_f(n_f,j)*v(n_f,j)*cp_f(n_f,j-1) - rho_f(n_f-1,j)*v(n_f-1,j)*cp_f(n_f-1,j-1))/dx...
            - (a_T(n_n,n_n-1) - max(0, -rho_f(n_f,j)*v(n_f,j)*cp_f(n_f,j-1)/dx));
        
        b_T(n_n) = Q_n(n_n,j) + (P(n_n,j)-P(n_n,j-1))/dt ...
                + v_n(n_n,j)*(P_f(n_f,j) - P_f(n_f-1,j))/dx ...
                + f(n_n)*rho(n_n,j)*abs(v_n(n_n,j))^3/(2*D) ...
                + rho(n_n,j-1)*cp(n_n,j-1)*T(n_n,j-1)/dt...
                + max(0, rho_f(n_f,j)*v(n_f,j)*cp_f(n_f,j-1)/dx)*T_f(n_f,j);
        
    
    
        % a_T(n_n) ??
    
        for i=2:n_n-1
            a_T(i,i-1) = -max( rho_f(i,j)*v(i,j)*cp_f(i,j-1)/dx, 0);
            a_T(i,i+1) = -max(   0, -rho_f(i+1,j)*v(i+1,j)*cp_f(i+1,j-1)/dx);
            a_T(i,i) = rho(i,j)*cp(i,j-1)/dt ...
                + (rho_f(i+1,j)*v(i+1,j)*cp_f(i+1,j-1) - rho_f(i,j)*v(i,j)*cp_f(i,j-1))/dx...
                - (a_T(i,i+1) + a_T(i,i-1));
            
            b_T(i) = Q_n(i,j) + (P(i,j)-P(i,j-1))/dt ...
                    + v_n(i,j)*(P_f(i+1,j) - P_f(i,j))/dx ...
                    + f(i)*rho(i,j)*abs(v_n(i,j))^3/(2*D) ...
                    + rho(i,j-1)*cp(i,j-1)*T(i,j-1)/dt;
        end
        
        T(:,j) = linsolve(a_T,b_T);
    
    end
    % THERMODYNAMIC PROPERTIES

    % Upwind scheme
    T_f(2:end-1,j) = (v(1:end-2,j) >= 0).*T(1:end-1,j) ...
                     + (v(1:end-2,j) <  0).*T(2:end,j);
    
    T_f(end,j) = T(end,j); % ASSUMING UPWIND SCHEME WITH FLOW FROM LEFT !!!!!!!!!!
                           % Implement other cases


    for i = 1:n_n  % PROPERTIES FROM P AND RHO          
        if strcmp(heat_transfer_model,'Adiabatic')
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),fluid);
            cp(i,j) = CP.PropsSI('C','P',P(i,j),'D',rho(i,j),fluid);
            h(i,j) = CP.PropsSI('H','P',P(i,j),'D',rho(i,j),fluid);     
            s(i,j) = CP.PropsSI('S','P',P(i,j),'D',rho(i,j),fluid);
            u(i,j) = CP.PropsSI('U','P',P(i,j),'D',rho(i,j),fluid);
            u_sonic_n(i,j) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),fluid);
            cp_f(i,j) = CP.PropsSI('C','P',P_f(i,j),'D',rho_f(i,j),fluid);
            h_f(i,j) = CP.PropsSI('H','P',P_f(i,j),'D',rho_f(i,j),fluid);
            s_f(i,j) = CP.PropsSI('S','P',P_f(i,j),'D',rho_f(i,j),fluid);
            u_sonic_f(i,j) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),fluid);
        elseif strcmp(heat_transfer_model,'Isothermal')
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),fluid);
            cp(i,j) = CP.PropsSI('C','T',T(i,j),'D',rho(i,j),fluid);
            h(i,j) = CP.PropsSI('H','T',T(i,j),'D',rho(i,j),fluid);     
            s(i,j) = CP.PropsSI('S','T',T(i,j),'D',rho(i,j),fluid);
            u(i,j) = CP.PropsSI('U','T',T(i,j),'D',rho(i,j),fluid);
            u_sonic_n(i,j) = CP.PropsSI('speed_of_sound','T',T(i,j),'D',rho(i,j),fluid);
            cp_f(i,j) = CP.PropsSI('C','T',T_f(i,j),'D',rho_f(i,j),fluid);
            h_f(i,j) = CP.PropsSI('H','T',T_f(i,j),'D',rho_f(i,j),fluid);
            s_f(i,j) = CP.PropsSI('S','T',T_f(i,j),'D',rho_f(i,j),fluid);
            u_sonic_f(i,j) = CP.PropsSI('speed_of_sound','T',T_f(i,j),'D',rho_f(i,j),fluid);
        elseif strcmp(heat_transfer_model,'Steady_state')
            error('Steady state under implementation.')
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),fluid);
            cp(i,j) = CP.PropsSI('C','T',T(i,j),'D',rho(i,j),fluid);
            h(i,j) = CP.PropsSI('H','T',T(i,j),'D',rho(i,j),fluid);     
            s(i,j) = CP.PropsSI('S','T',T(i,j),'D',rho(i,j),fluid);
            u(i,j) = CP.PropsSI('U','T',T(i,j),'D',rho(i,j),fluid);
            u_sonic_n(i,j) = CP.PropsSI('speed_of_sound','T',T(i,j),'D',rho(i,j),fluid);
            cp_f(i,j) = CP.PropsSI('C','T',T_f(i,j),'D',rho_f(i,j),fluid);
            h_f(i,j) = CP.PropsSI('H','T',T_f(i,j),'D',rho_f(i,j),fluid);
            s_f(i,j) = CP.PropsSI('S','T',T_f(i,j),'D',rho_f(i,j),fluid);
            u_sonic_f(i,j) = CP.PropsSI('speed_of_sound','T',T_f(i,j),'D',rho_f(i,j),fluid);
        elseif strcmp(heat_transfer_model,'Kiuchi')
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),fluid);
            cp(i,j) = CP.PropsSI('C','T',T(i,j),'D',rho(i,j),fluid);
            h(i,j) = CP.PropsSI('H','T',T(i,j),'D',rho(i,j),fluid);     
            s(i,j) = CP.PropsSI('S','T',T(i,j),'D',rho(i,j),fluid);
            u(i,j) = CP.PropsSI('U','T',T(i,j),'D',rho(i,j),fluid);
            u_sonic_n(i,j) = CP.PropsSI('speed_of_sound','T',T(i,j),'D',rho(i,j),fluid);
            cp_f(i,j) = CP.PropsSI('C','T',T_f(i,j),'D',rho_f(i,j),fluid);
            h_f(i,j) = CP.PropsSI('H','T',T_f(i,j),'D',rho_f(i,j),fluid);
            s_f(i,j) = CP.PropsSI('S','T',T_f(i,j),'D',rho_f(i,j),fluid);
            u_sonic_f(i,j) = CP.PropsSI('speed_of_sound','T',T_f(i,j),'D',rho_f(i,j),fluid);
        else
            error('Heat transfer model not identified for thermodynamic properties calculations')
        end
    end
    
    cp_f(n_f,j) = CP.PropsSI('C','P',P_f(n_f,j),'D',rho_f(n_f,j),fluid);
    h_f(n_f,j) = CP.PropsSI('H','P',P_f(n_f,j),'D',rho_f(n_f,j),fluid);
    s_f(n_f,j) = CP.PropsSI('S','P',P_f(n_f,j),'D',rho_f(n_f,j),fluid);
    u_sonic_f(n_f,j) = CP.PropsSI('speed_of_sound','P',P_f(n_f,j),'D',rho_f(n_f,j),fluid);
    

    % CAES process
    % Identify and update operational state (Charging, discharging, idle)
    if strcmp(Process,'Charging_L')
        % if strcmp(L_bound,'Inlet') & P(ceil(n_n/2),j) >= P_max
        if strcmp(L_bound,'Inlet') & P(1,j) >= P_max
            % Charging from L boundary
            % v(1,j+1:end) = 0;
            L_bound = 'Wall';
            t_shut_off = (j-1)*dt; 
            tol = 1e-8;
        end
    elseif strcmp(Process,'Discharging_L')
        if strcmp(L_bound,'M_const') & P(1,j) <= P_min
        % if strcmp(L_bound,'M_const') & P(ceil(n_n/2),j) <= P_min
            % Discharging from L boundary
            L_bound = 'Wall';
            t_shut_off = (j-1)*dt;
            tol = 1e-8;
        end
    elseif strcmp(Process,'Charging_R')
        % if strcmp(R_bound,'Inlet') & P(end-1,j) >= P_max
        if strcmp(R_bound,'Inlet') & P(ceil(n_n/2),j) >= P_max
            % Charging from R boundary
            % v(1,j+1:end) = 0;
            R_bound = 'Wall';
            t_shut_off = (j-1)*dt;
        end
    elseif strcmp(Process,'Discharging_R')
        % if strcmp(R_bound,'M_const') & P(end,j) <= P_min
        if strcmp(R_bound,'M_const') & P(ceil(n_n/2),j) <= P_min
            % Discharging from R boundary
            R_bound = 'Wall';
            t_shut_off = (j-1)*dt;
        end
    elseif strcmp(Process,'Cycle_L')
        if strcmp(stage,'Charging') & P(1,j) >= P_max
            % Charging from L boundary
            % v(1,j+1:end) = 0;
            L_bound = 'Wall';
            stage = 'idle_charg';
            % t_ch = (j-1)*dt;
        % elseif strcmp(stage,'idle_charg') & (std(P(:,j)./P(:,j)) <= 0.1 | t(j) >= Dt_charg)
            j_charg_end = j;
            disp(strcat('Charging completed, ',string(timeofday(datetime) ) ))
        elseif strcmp(stage,'idle_charg') & t(j) >= Dt_charg
            L_bound = 'M_const';
            j_disch = j; % index of start of discharge
            m_L = -100;
            stage = 'Discharging';
            disp(strcat('Discharging started, ',string(timeofday(datetime)) ) )
        elseif strcmp(stage,'Discharging') & P(1,j) <= P_min
            L_bound = 'Wall';
            stage = 'idle_disch';
            j_disch_end = j;
            disp(strcat('Discharging completed, ',string(timeofday(datetime)) ) )
        end
        stage_hist = [stage_hist;stage];
    elseif strcmp(Process,'Cycle_R')
    elseif strcmp(Process,'NCycle_L')
    elseif strcmp(Process,'NCycle_R')
    elseif strcmp(Process,'Kiuchi')


        % warning('Kiuchi boundary update procedure under implementation.')
        
        
        if strcmp(stage,'Idle_Kiuchi') & (t(j-1) < 10*60 & t(j) >= 10*60)
            
            L_bound = 'P_const';

            R_bound = 'Vol_const';

            stage = 'Active_Kiuchi';
        elseif strcmp(stage,'Active_Kiuchi') & (t(j-1) < 30*60 & t(j) >= 30*60)
            L_bound = 'Wall';
            R_bound = 'Wall';
            stage = 'Idle';
        end
        stage_hist = [stage_hist;stage];
    else
        error('Process not identified.')
    end


    bound_hist = [bound_hist; string(L_bound) string(R_bound)];
    
    m(j) = sum(rho(:,j)*A_h*dx);
    % E(j) = sum(rho(:,j)*A_h*dx.*u(:,j));
    E(j) = sum(rho(:,j)*A_h*dx.*( u(:,j) + v_n(:,j).^2/2 ) );
    
    m_n(:,j) = rho(:,j)*A_h*dx;
    % E_n(:,j) = rho(:,j)*A_h*dx.*u(:,j);
    E_n(:,j) = rho(:,j)*A_h*dx.*( u(:,j) + v_n(:,j).^2/2 );
    
    % if rem(t(j),1)==0
    %     disp(strcat('t = ',num2str(t(j)),'s'))
    % end
    if t(j-1) < Dt/5 && t(j) >= Dt/5
        disp(strcat('20% of time steps.',string(timeofday(datetime) ) ) )
    elseif t(j-1) < 2*Dt/5 && t(j) >= 2*Dt/5
        disp(strcat('40% of time steps.',string(timeofday(datetime))))
    elseif t(j-1) < 3*Dt/5 && t(j) >= 3*Dt/5
        disp(strcat('60% of time steps.',string(timeofday(datetime))))
    elseif t(j-1) < 4*Dt/5 && t(j) >= 4*Dt/5
        disp(strcat('80% of time steps.',string(timeofday(datetime) )))
    end

end

elapsedTime = toc;
if elapsedTime/60 > 5
    disp(strcat('Elapsed time: ',num2str(floor(elapsedTime/60)),' min'));
else
    disp(strcat('Elapsed time: ',num2str(floor(elapsedTime)),' s'));
end

dm = zeros(n_t,1);
dE = zeros(n_t,1);

% Mass, energy and exergy increments at each time step
if strcmp(Process,'Charging_L')
    % dm = rho_L*v(1,:)'*A_h*dt;
    % dE = rho_L*v(1,:)'*A_h*cp_L*T_L*dt;
    % dX = rho_L*v(1,:)'*A_h*(h_L - h_o - T_o*(s_L - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible

    % Sliding pressure test 
    dm = rho_f(1,:)'.*v(1,:)'*A_h*dt;
    % dE = rho_f(1,:)'.*v(1,:)'.*A_h.*h_f(1,:)'*dt;
    dE = rho_f(1,:)'.*v(1,:)'.*A_h.*( h_f(1,:)' + v(1,:)'.^2/2 )*dt;
    dX = rho_f(1,:)'.*v(1,:)'*A_h.*(h_f(1,:)' - h_o - T_o*(s_f(1,:)' - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible

elseif strcmp(Process,'Discharging_L')
    dm = rho_f(1,:)'.*v(1,:)'*A_h*dt;
    dE = rho_f(1,:)'.*v(1,:)'.*A_h.*(h_f(1,:)' + v(1,:)'.^2/2 )*dt;
    dX = rho_f(1,:)'.*v(1,:)'*A_h.*(h_f(1,:)' - h_o - T_o*(s_f(1,:)' - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible

    dX_m = rho_f(1,:)'.*v(1,:)'*A_h;
    dX_h = h_f(1,:)' - h_o;
    dX_s = T_o*(s_f(1,:)' - s_o);
elseif strcmp(Process,'Charging_R')
    dm = -rho_R*v(end,:)'*A_h*dt;
    dE = -rho_R*v(end,:)'*A_h*cp_R*T_R*dt;
    dX = -rho_R*v(end,:)'*A_h*(h_R - h_o - T_o*(s_R - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible

elseif strcmp(Process,'Discharging_R')
    dm = -rho_f(end,:)'.*v(end,:)'*A_h*dt;
    dE = -rho_f(end,:)'.*v(end,:)'*A_h.*cp_f(end,:)'.*T_f(end,:)'*dt;
    dX = -rho_f(end,:)'.*v(end,:)'*A_h.*(h(:,end)' - h_o - T_o*(s(:,end)' - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible

elseif strcmp(Process,'Cycle_L')
    dm = rho_f(1,:)'.*v(1,:)'*A_h*dt;
    % dE = rho_f(1,:)'.*v(1,:)'.*A_h.*h_f(1,:)'*dt;
    dE = rho_f(1,:)'.*v(1,:)'.*A_h.*( h_f(1,:)' + v(1,:)'.^2/2 )*dt;
    dX = rho_f(1,:)'.*v(1,:)'*A_h.*(h_f(1,:)' - h_o - T_o*(s_f(1,:)' - s_o))*dt;

elseif strcmp(Process,'Cycle_R')
    error('Cycle_R process not implemented yet ')
elseif strcmp(Process,'NCycles_L')
    error('NCycle_L process not implemented yet ')
elseif strcmp(Process,'NCycles_R')
    error('NCycle_R process not implemented yet ')
elseif strcmp(Process,'Idle')
    dm = 0;
    dE = Q;
    dX = 0;
elseif strcmp(Process,'Kiuchi')
    % Not required for Kiuchi
    dm = 0;
    dE = 0;
    dX = 0;
    % error('Kiuchi energy, mass and exergy variation not implemented yet.')
else
    error('Process not identified.')
end

% Energy and mass balance: Property at t=0 + cummulative sum of mass/energy
% that enters the control volume
m_bal = m(1) + cumsum(dm);
E_bal = E(1) + cumsum(dE);

% Exergy (Assuming compressed air temperature equal to environmental temperature
% X = P*(A_h*L).*(P_amb./P - 1 + log(P./P_amb))./(1e6*3600);            % Pipeline Exergy [MWh] 
% X_min = P_0*(A_h*L).*(P_amb./P_0 - 1 + log(P_0./P_amb))/(1e6*3600);   % Exergy when discharged [MWh] 
% X = P*(A_h*dx).*(P_amb./P - 1 + log(P./P_amb))./(1e6*3600);           % Pipeline Exergy [MWh]
% X_min = P_0*(A_h*dx).*(P_amb./P_0 - 1 + log(P_0./P_amb))/(1e6*3600);  % Exergy when discharged [MWh] 

% Total Exergy in the pipeline at each time step [MWh]
X = sum(m_n.*( u - u_o + P_o*R*(T./P - T_o/P_o) - T_o*(s - s_o) ) )./(1e6*3600);
% Exergy at t = 0s [MWh] 
X_0 = m(1)*(u_0 - u_o + P_o*R*(T_0./P_0 - T_o/P_o) - T_o*(s_0 - s_o) )/(1e6*3600); 


X_net = X - X_0;                                                 % Exergy between current state and discharged state (assuming whole pipeline at P_min)
DeltaX_flow = sum(dX)/(1e6*3600);
DeltaX_st = sum(X_net(:,end));
if strcmp(Process,'Charging_L')
    etaX_Ch = DeltaX_st/DeltaX_flow
elseif strcmp(Process,'Discharging_L')
    etaX_Disch = DeltaX_flow/DeltaX_st
end

% figure('Color',[1 1 1])
% if strcmp(Process,'Cycle_L') | strcmp(Process,'Cycle_R')
%     X_0_disch = sum(m_n(:,j_disch).*(u(:,j_disch) - u_o + P_o*R*(T(:,j_disch)./P(:,j_disch) - T_o/P_o) - T_o*(s(:,j_disch) - s_o) ))/(1e6*3600);  % Exergy at t = 0s [MWh] 
%     if Figures
%         plot(t,X,t,[X_0+cumsum(dX(t<Dt_charg))./(1e6*3600);X_0_disch+cumsum(dX(t>=Dt_charg))./(1e6*3600)])
%         legend('X','X_0 + dX')
%     end
% else
%     if Figures
%         plot(t,X,t,X_0+cumsum(dX)./(1e6*3600))
%         legend('X','X_0 + dX')
%     end
% end


if strcmp(Process,'Cycle_L') | strcmp(Process,'Cycle_R')
    m_bal(t>=Dt_charg) = m(j_disch) + cumsum(dm(t>=Dt_charg));
    E_bal(t>=Dt_charg) = E(j_disch) + cumsum(dE(t>=Dt_charg));
    
    if Figures
        figure('color',[1 1 1]);plot(t,m)
        hold on; plot(t,m_bal)
        legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
        figure('color',[1 1 1]);plot(t,E./(1e6*3600))
        hold on; plot(t,E_bal(1:end)./(1e6*3600))
        grid on;
        legend('E [MWh]','$E_o + \dot{m} \Delta E [MWh]$','Interpreter','latex')
    
        figure('color',[1 1 1]);
        subplot(1,2,1)
        t_charg = t(t<Dt_charg);
        plot(t_charg,X(t<Dt_charg),t_charg,X(1)+cumsum(dX(t<Dt_charg))./(1e6*3600))
        title('Charging')
        ylim([2.9 5])
        subplot(1,2,2)
        t_disch = t(t>=Dt_charg) - t_charg(end);
        plot(t_disch,X(t>=Dt_charg),t_disch,X(j_disch)+cumsum(dX(t>=Dt_charg))./(1e6*3600))
        title('Discharging')
        ylim([2.9 5])
    
        % figure('color',[1 1 1]);
        % subplot(1,2,1)
        % title('Charging')
        % plot(t(t<Dt_charg),X(t<Dt_charg),t(t<Dt_charg),X(1)+cumsum(dX(t<Dt_charg))./(1e6*3600))
        % subplot(1,2,2)
        % title('Discharging')
        % plot(t(t>=Dt_charg),X(t>=Dt_charg),t(t>=Dt_charg),X(j_disch)+cumsum(dX(t>=Dt_charg))./(1e6*3600))
    end
elseif strcmp(Process,'Charging_L')
    
    if Figures
        figure('color',[1 1 1]);plot(t,m - m(1))
        hold on; plot(t,m_bal - m(1))
        grid on
        legend('$\Delta m_{storage}$','$m_{in}$','Interpreter','latex')

        figure('color',[1 1 1]);plot(t,(E - E(1))./(1e6*3600))
        hold on; plot(t,(E_bal(1:end) - E(1))./(1e6*3600))
        grid on 
        legend('$\Delta E_{storage}$ [MWh]','$E_{in}$ [MWh]','Interpreter','latex')
    end
elseif strcmp(Process,'Discharging_L')
    if Figures
        figure('color',[1 1 1]);plot(t,m(1)-m)
        hold on; plot(t,m(1) - m_bal)
        grid on
        legend('$\Delta m_{storage}$','$m_{out}$','Interpreter','latex')
        % legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
        
        figure('color',[1 1 1]);plot(t,(E(1) - E)./(1e6*3600))
        hold on; plot(t,(E(1) - E_bal(1:end))./(1e6*3600))
        grid on 
        legend('$\Delta E_{storage}$ [MWh]','$E_{out}$ [MWh]','Interpreter','latex')
    end
else
    if Figures
        figure('color',[1 1 1]);plot(t,m)
        hold on; plot(t,m_bal)
        legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
        figure('color',[1 1 1]);plot(t,E./(1e6*3600))
        hold on; plot(t,E_bal(1:end)./(1e6*3600))
        grid on 
        legend('E [MWh]','$E_o + \dot{m} \Delta E [MWh]$','Interpreter','latex')
    end
end


% Analysing velocity profiles at first node (at each face and at the centre
% of the node)
% figure('Color',[1 1 1])
% plot(t,v(1,:))
% hold on
% plot(t,v_n(1,:))
% plot(t,v(2,:))
% legend('v(1)','v_n(1)','v(2)')
% xlabel('t')
% ylabel('v')

if Figures
    figure('Color',[1 1 1]);
    plot(abs(mean(error_hist)))
    xlabel('Iteration')
    ylabel('Absolute mean of the error')
    grid on
end

m2 = sum(rho(2:end,:)*A_h*dx);
E2 = sum(rho(2:end,:)*A_h*dx.*( u(2:end,:) + v_n(2:end,:).^2/2 ) );

dm2 = rho_f(2,:)'.*v(2,:)'*A_h*dt;
    % dE = rho_f(1,:)'.*v(1,:)'.*A_h.*h_f(1,:)'*dt;
dE2 = rho_f(2,:)'.*v(2,:)'.*A_h.*( h_f(2,:)' + v(2,:)'.^2/2 )*dt;

m2_bal = m2(1) + cumsum(dm2);
E2_bal = E2(1) + cumsum(dE2);

if Figures
    figure;
    yyaxis left
    plot(t,m2,t,m2_bal)
    ylabel('mass [kg]')
    yyaxis right
    plot(t,E2./(1e6*3600),t,E2_bal./(1e6*3600))
    ylabel('Energy [MWh]')
    legend('m','m_{bal}','E','E_{bal}')
    
    if strcmp(Process,'Charging_L')
        figure('color',[1 1 1]);plot(t./3600,(E2 - E2(1))./(1e6*3600))
        hold on; plot(t./3600,(E2_bal(1:end) - E2(1))./(1e6*3600))
        grid on 
        legend('$\Delta E_{storage}$ [MWh]','$E_{in}$ [MWh]','Interpreter','latex')
    elseif strcmp(Process,'Discharging_L')
        figure('color',[1 1 1]);plot(t./3600,(E2(1) - E2)./(1e6*3600))
        hold on; plot(t./3600,(E2(1) - E2_bal(1:end))./(1e6*3600))
        grid on 
        legend('$\Delta E_{storage}$ [MWh]','$E_{out}$ [MWh]','Interpreter','latex')
    end

    figure
    plot(t,(E2 - E2_bal')./E2)
    xlabel(' t [s]')
    ylabel('Relative difference between E and E_{bal}')
end

XX = sum(m_n(2:end,:).*( u(2:end,:) - u_o + P_o*R*(T(2:end,:)./P(2:end,:) - T_o/P_o) - T_o*(s(2:end,:) - s_o) ) )./(1e6*3600);

dXX = rho_f(2,:)'.*v(2,:)'*A_h.*(h_f(2,:)' - h_o - T_o*(s_f(2,:)' - s_o))*dt./(1e6*3600); % [MWh]

XX_bal = XX(1) + cumsum(dXX);

if strcmp(Process,'Charging_L')
    etaX_stor = ( XX(end) - XX(1) ) / sum(dXX)
elseif strcmp(Process,'Discharging_L')
    etaX_stor = sum(dXX) / ( XX(end) - XX(1) )
elseif strcmp(Process, 'Cycle_L')
    error('eta_X calculation not implemented for cycles yet. ')
    
elseif strcmp(Process,'Kiuchi')
    etaX_stor = 0;
    figure('Color',[1 1 1])
    plot(t./60,rho_f(1,:).*v(1,:)*A_h*3600/rho_st)
    hold  on
    plot(t./60,rho_f(end,:).*v(end,:)*A_h*3600/rho_st)
    legend('Face 1','Face n')
    ylim([-1e5 4e5])
    ylabel('Volumetric flow rate [scmh]')

    figure('Color',[1 1 1])
    plot(t./60,rho(1,:).*v_n(1,:)*A_h*3600/rho_st,'b')
    hold  on
    plot(t./60,rho(end,:).*v_n(end,:)*A_h*3600/rho_st,'r')
    legend('Node 1','Node n')
    ylabel('Volumetric flow rate [scmh]')

    figure('Color',[1 1 1])
    plot(t./60,rho_f(1,:).*v(1,:)*A_h)
    hold  on
    plot(t./60,rho_f(end,:).*v(end,:)*A_h)
    legend('Face 1','Face n')
    ylabel('Mass flow rate [kg/s]')
else
    error('Process not found or eta not implemented for the process.')
end

if Save_data
    if strcmp(Process,'Charging_L')
        filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'-',num2str(P_min/1e6),'MPa_m_',num2str(m_in));
        % filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'-',num2str(P_min/1e6),'MPa_m_',num2str(m_in),'_',string(year(time)),'_',string(month(time)),'_',string(day(time)));
    elseif strcmp(Process,'Discharging_L')
        filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'-',num2str(P_min/1e6),'MPa_m_',num2str(m_out));
        % filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'-',num2str(P_min/1e6),'MPa_m_',num2str(m_out),'_',string(year(time)),'_',string(month(time)),'_',string(day(time)));
    elseif strcmp(Process,'Cycle_L')
        filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'-',num2str(P_min/1e6),'MPa_m_',num2str(m_in));
    elseif strcmp(Process,'Kiuchi')
        time = now;
        filename = strcat(Process,'_',string(month(time)),'_',string(day(time)));
    else
        error(strcat('Save data for process ',Process,' not implemented'))
    end
    % filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'MPa_m_',num2str(m_in),'_',string(year(time)),'_',string(month(time)),'_',string(day(time)));
    save(filename)
end

disp(strcat('Boundary condition: ',L_bound))
%% Previous tests
if Figures
    % figure('color',[1 1 1]);plot(t,(m_bal' - m)./m)
    % title('Difference in mass')
    % figure('color',[1 1 1]);plot(t,(E_bal' - E)./E)
    % title('Difference in energy')
    
    % Maximum temperature point variation
    x_maxT = zeros(length(t),1);
    Tmax = zeros(length(t),1);
    for i=1:length(t)
        [Tmax(i), ind_max] = max(T(:,i));
        x_maxT(i) = x(ind_max);
    end
    figure('Color',[1 1 1])
    yyaxis left
    plot(t,x_maxT./1e3)
    ylabel('Position [km]')
    yyaxis right
    plot(t,Tmax);
    ylabel('T_{max} [K]')
end
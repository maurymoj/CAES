clear 
% close all
% clc
datetime
% profile on
tic
CP = py.importlib.import_module('CoolProp.CoolProp'); % Simplifies coolprop calls

%----------------------- PROBLEM PARAMETERS -------------------------%
% Pipeline parameters
% Parameters from Kiuchi (1994)
% L = 5000;
% D = ;
% L = 70000; % Reference case 70 km length, 0.9 m diam
% D = 0.9;
% L = 213333; % Length for eq. vol with Huntorf
L = 100000;
D = 0.9;

% Cavern parameters
% L = 300; % Volume shape similar to Huntorf
% D = 24;
% L = 35; % Roughly same volume as 70 km, 0.9 D pipeline
% D = 40;

% Total simulation time
% Dt = 360; % Charging time for 5 km pipeline
% Dt = 720; % Charging + discharging time (5km 0.5m pipeline)
% Dt_charg = 360;
% Dt = 1200; % Charging L=10 km, d=0.9 m pipeline (360s = 3.44 MWh,1054 elapsed time)
Dt = 4*3600;

eps = 0.04e-3; % Absolute roughness 0.04 mm

% Ambient conditions
P_amb = 101325;
T_amb = 273.15 + 25;

T_ground = 278;

% System operational limits
P_max = 9e6;
% P_min = 4.3e6; % Huntorf
% P_min = 4e6;
P_min = P_max - 3e6;

% Charging rate
m_in = 200;
T_in = 273.15+60;

% Discharging rate
m_out = 200;

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
Process = 'Charging_L';

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
    T_0 = 273.15 + 25;
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

    L_bound = 'Inlet'; % OR M_CONST ??
    R_bound = 'Wall';

    stage = 'Charging';
    stage_hist = stage;
elseif strcmp(Process,'Cycle_R')
elseif strcmp(Process,'NCycles_L')
elseif strcmp(Process,'NCycles_R')
elseif strcmp(Process,'Idle')
    P_0 = P_min;
    T_0 = T_amb;
    v_0 = 0;

    L_bound = 'Wall';
    R_bound = 'Wall';
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
    % m_A = -417; % Huntorf, sign indicates flow direction
    m_L = m_out;
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
    % m_B = 417; % Huntorf    
    m_R = -m_out;
end                     

% heat_transfer_model 
% 'Adiabatic'
% ' Steady_state'
heat_transfer_model = 'Adiabatic';

k_pipe = 45.3; % [W/m K] - Wen et al. (2023) "Heat Transfer Model of Natural Gas Pipeline Based on Data Feature Extraction and First Principle Models"
k_int_coating = 0.52;
k_ext_coating = 0.4;
cp_pipe = 500; % [J/kg K]

thickness_pipe = 15.9e-3; % Still based on Wen, but in accordance with the regulations in the UK
thickness_int_coating = 0.5e-3;
thickness_ext_coating = 3e-3;

R = 287;
g = 9.81;
theta = 0;

Figures = 0;

%--------------------- SIMULATION PARAMETERS ------------------------%

if strcmp(simType,'CAESPipe')
    n_nodes = 40;
    max_iter = 20;
    % dx = L/(100-1); % CAESPipe n_nodes = 100 and dt = 0.25 provided low
    % residual
    dx = L/(n_nodes-1); 
    dt_max = dx/400; % 400 is representative of the sound speed, 
                     % it is higher than the maximum sound speed reached in the pipeline 
                     % to achieve a conservative value
    dt = (Dt/ceil(Dt/dt_max)) % division of Dt in an integer number of intervals
                            % with dt smaller than dt_max
    if dx/400 < dt % 400 upper limit for the speed of sound
        warning('dt > time needed for pressure wave to cross a node')
    end
    
    if strcmp(Process,'Discharging_R') || strcmp(Process,'Discharging_L')
        tol = 1e-5; % CAESPipe discharging
    elseif strcmp(Process,'Charging_R') || strcmp(Process,'Charging_L')
        tol = 1e-6; % CAESPipe Charging
    elseif strcmp(Process,'Cycle_L')
        tol = 1e-5;
    elseif strcmp(Process,'Idle')
        tol = 1e-5;
    else
        error('Process not found.')
    end

elseif strcmp(simType,'CAESCav')
    n_nodes = 5;
    max_iter = 20;
    dx = L/(n_nodes-1);
    dt_max = dx/400;
    dt = (Dt/ceil(Dt/dt_max))
    tol = 1e-7;
else
    warning('Cant identify simulation type.')
end

% Friction calculation equation
% 'Nikuradse'
% 'Colebrook'; % No significant difference and adds calculation time
friction_model = 'Nikuradse';

% Tuning

% Under-relaxation (1 means no under-relaxation)
alpha = 0.5;
alpha_P = alpha ;  % Pressure under-relaxation factor
alpha_v = alpha ;  % velocity under-relaxation factor
alpha_rho = alpha ;  % Density under-relaxation factor
%---------------------- ARRAYS INITIALIZATION ----------------------%
t = 0:dt:Dt;

n_n = L/dx;       % number of nodes
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

% Node initialization
v_n = zeros(n_n,n_t);
P = zeros(n_n,n_t);
T = zeros(n_n,n_t);
rho = zeros(n_n,n_t);
h = zeros(n_n,n_t);
s = zeros(n_n,n_t);
u = zeros(n_n,n_t);
cp = zeros(n_n,n_t);
% equation and compare to using only mass and momentum eqs
Q = zeros(1,n_t);
Q_n = zeros(n_n,n_t);

% Sanity check
m = zeros(1,n_t); % Total mass in pipeline
E = zeros(1,n_t); % Total energy in pipeline

m_n = zeros(n_n,n_t); % Mass at each node
E_n = zeros(n_n,n_t); % Energy at each node

%------------------ SETTING INITIAL CONDITIONS ---------------------%
% Basic calculations
epsD = eps/D;   % Relative roughness
A_h = pi*D^2/4; % Cross-section area

% Dead state parameters
P_o = P_amb;
T_o = T_amb;
h_o = CP.PropsSI('H','P',P_o,'T',T_o,'Air');
u_o = CP.PropsSI('U','P',P_o,'T',T_o,'Air');
s_o = CP.PropsSI('S','P',P_o,'T',T_o,'Air');

% Initial conditions at nodes (i -> x, j -> t)
u_0 = CP.PropsSI('U','P',P_0,'T',T_0,'Air');
s_0 = CP.PropsSI('S','P',P_0,'T',T_0,'Air');

% Thermodynamic properties over all solution space
P(:,1) = P_0*ones(n_n,1);
T(:,1) = T_0*ones(n_n,1);
rho(:,1) = CP.PropsSI('D','P',P(1,1),'T',T(1,1),'Air');
cp(:,1) = CP.PropsSI('C','P',P(1,1),'T',T(1,1),'Air'); 
h(:,1) = CP.PropsSI('H','P',P(1,1),'T',T(1,1),'Air'); 
s(:,1) = CP.PropsSI('S','P',P(1,1),'T',T(1,1),'Air'); 
u(:,1) = CP.PropsSI('U','P',P(1,1),'T',T(1,1),'Air'); 

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
    rho_L = CP.PropsSI('D','P',P_L,'T',T_L,'Air');
    cp_L = CP.PropsSI('C','P',P_L,'T',T_L,'Air');
    h_L = CP.PropsSI('H','P',P_L,'T',T_L,'Air');
    s_L = CP.PropsSI('S','P',P_L,'T',T_L,'Air');

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
end


% R boundary conditions
if strcmp(R_bound,'Inlet')
    rho_R = CP.PropsSI('D','P',P_R,'T',T_R,'Air');
    cp_R = CP.PropsSI('C','P',P_R,'T',T_R,'Air');
    h_R = CP.PropsSI('H','P',P_R,'T',T_R,'Air');
    s_R = CP.PropsSI('S','P',P_R,'T',T_R,'Air');

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
cp_f(1,1) = CP.PropsSI('C','P',P_f(1,1),'D',rho_f(1,1),'Air');
h_f(1,1) = CP.PropsSI('H','P',P_f(1,1),'D',rho_f(1,1),'Air');
s_f(1,1) = CP.PropsSI('S','P',P_f(1,1),'D',rho_f(1,1),'Air');

cp_f(end,1) = CP.PropsSI('C','P',P_f(end,1),'D',rho_f(end,1),'Air');
h_f(end,1) = CP.PropsSI('H','P',P_f(end,1),'D',rho_f(end,1),'Air');
s_f(end,1) = CP.PropsSI('S','P',P_f(end,1),'D',rho_f(end,1),'Air');

for i = 1:n_n  % PROPERTIES FROM P AND RHO          
    % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),'Air');
    cp(i,1) = CP.PropsSI('C','P',P(i,1),'D',rho(i,1),'Air');
    h(i,1) = CP.PropsSI('H','P',P(i,1),'D',rho(i,1),'Air');
    s(i,1) = CP.PropsSI('S','P',P(i,1),'D',rho(i,1),'Air');
    u(i,1) = CP.PropsSI('U','P',P(i,1),'D',rho(i,1),'Air');
    cp_f(i,1) = CP.PropsSI('C','P',P_f(i,1),'D',rho_f(i,1),'Air');
    h_f(i,1) = CP.PropsSI('H','P',P_f(i,1),'D',rho_f(i,1),'Air');
    s_f(i,1) = CP.PropsSI('S','P',P_f(i,1),'D',rho_f(i,1),'Air');
end

m(1) = sum(rho(:,1)*A_h*dx);
E(1) = sum(rho(:,1)*A_h*dx.*u(:,1));

m_n(:,1) = rho(:,1)*A_h*dx;
E_n(:,1) = rho(:,1)*A_h*dx.*u(:,1);

count = zeros(n_t,1);

% Friction factor based on Nikuradse
f_guess = (2*log10(1/epsD)+1.14)^(-2);

% resFig = figure;
error_hist = [];
bound_hist = [string(L_bound) string(R_bound)];

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
    count(j) = 0;
    error_P = 10;

    % Momentum and mass balance loop
    while count(j) < max_iter && max(abs(error_P)) > tol
        % Under-relaxed corrections
        P(:,j) = P(:,j) + alpha_P*P_corr;
        rho(:,j) = alpha_rho*(rho(:,j) + rho_corr) ...
            + (1-alpha_rho)*rho(:,j);

        v(:,j) = alpha_v*(v_star + v_corr) ...
            + (1-alpha_v)*v_star;
        % v(1:end-1,j) = alpha_v*(v_star(1:end-1) + v_corr(1:end-1)) ...
        % + (1-alpha_v)*v_star(1:end-1);

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
            rho_f(1,j) = CP.PropsSI('D','P',P_f(1,j),'T',T_f(1,j),'Air');
            v(1,j) = m_L/(rho_f(1,j)*A_h); % Sliding pressure
            
            v_n(1,j) = m_L/(rho(1,j)*A_h);

        elseif strcmp(L_bound,'Wall')
            v(1,j) = 0;
            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);

            v_n(1,j) = (v(1,j) >= 0).*v(1,j) ...
            +      (v(1,j) <  0).*v(2,j);
        elseif strcmp(L_bound,'P_const')
            % P(1,:) = P_L;
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

            % P_f(1,j) = 2*P_f(2,j) - P_f(3,j);
            % T_f(1,j) = 2*T_f(2,j) - T_f(3,j);
            % rho_f(1,j) = 2*rho_f(2,j) - rho_f(3,j);
 
            v(1,j) = m_L/(rho_f(1,j)*A_h);

            % v_n(1,j) = m_L/(rho(1,j)*A_h);
            % Zero-gradient assumption
            % v_n(1,j) = v(1,j);
            % Upwind scheme
            v_n(1,j) = (v(1,j) >= 0).*v(1,j) ...
                +      (v(1,j) <  0).*v(2,j);
            % Central scheme
            % v_n(1,j) = (v(1,j) + v(2,j))/2;
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
            rho_f(end,j) = CP.PropsSI('D','P',P_f(end,j),'T',T_f(end,j),'Air');
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
        
        % Properties in nodes to faces
        P_f(2:end-1,j) = (v(2:end-1,j) >= 0).*P(1:end-1,j) ...
            +            (v(2:end-1,j) <  0).*P(2:end,j);
        rho_f(2:end-1,j) = (v(2:end-1,j) >= 0).*rho(1:end-1,j) ...
            +            (v(2:end-1,j) <  0).*rho(2:end,j);


        % Properties
        u_sonic_f = zeros(n_f,1);
        drho_dP_f = zeros(n_f,1);
        nu = zeros(n_f,1);

        u_sonic_n = zeros(n_n,1);
        drho_dP_n = zeros(n_n,1);

        for i = 1:n_n  % PROPERTIES FROM P AND RHO
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),'Air');

            u_sonic_f(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),'Air');
            if strcmp(friction_model,'Colebrook')
                nu = CP.PropsSI('viscosity','P',P_f(i,j),'D',rho_f(i,j),'Air');
            end

            drho_dP_f(i) = 1/(u_sonic_f(i)^2);

            u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),'Air');
            drho_dP_n(i) = 1/(u_sonic_n(i)^2);
        end

        u_sonic_f(end) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'D',rho_f(end,j),'Air');
        drho_dP_f(end) = 1/(u_sonic_f(end)^2);
        if strcmp(friction_model,'Colebrook')
            nu(end) = CP.PropsSI('viscosity','P',P_f(end,j),'D',rho_f(end,j),'Air');
        end


        %---------- v* calculation - Momentum control volume ---------------------%
        
        % Friction factor based on Nikuradse - IMPLEMENT COLEBROOK EQUATION
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
        end
        
        

        % Matrix/vector initialization (only need to solve for the inner faces -
        % boundary faces are dealt with the boundary conditions)
        a_M = zeros(n_f);
        b_M = zeros(n_f,1);
        d_M = zeros(n_f,1);
        B_M = zeros(n_f,1);
        
        a_M(1,1) = 1; % Value for v at the first momentum volume is 
                    % known from boundary conditions
        B_M(1) = v(1,j);
        
        a_M(n_f,n_f) = 1;
        B_M(n_f) = v(n_f,j);
        

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
        
        if strcmp(L_bound,'M_const') %|| strcmp(L_bound,'Inlet')
            % P_corr(1) = 0;
        elseif strcmp(R_bound,'M_const')
            P_corr(end) = 0;
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

        % v_corr(end)       =  % Velocity corrected based on mass balance (m_in =
        % m_out)
        error_P = (P_corr./P(:,j));
        error_hist = [error_hist error_P];
        
        % if rem(count(j),2) == 0 & Figures
        %     figure(resFig)
        %     plot(mean(abs(error_hist)))
        %     drawnow limitrate
        % end

        % error_rho = (rho_corr./rho(:,j));
        % error_v = (v_corr./v(:,j));
        count(j) = count(j)+1;  
    end

    % ENERGY BALANCE
    % ASSUMING v>0
    a_T = zeros(n_n,n_n);
    b_T = zeros(n_n,1);

    if strcmp(heat_transfer_model,'Adiabatic')
        Q(j) = 0;
    elseif strcmp(heat_transfer_model,'Sukhov')
    elseif strcmp(heat_transfer_model,'Steady_state')
        % neglecting effect of coating
        D_ext = D + thickness_pipe; 
        U = dx*k_pipe/(D*log(D_ext/D));
        Q_n(:,j) = -4*U/D*(T(:,j-1) - T_ground);
        Q(j) = sum(Q_n(:,j));
    elseif strcmp(heat_transfer_model,'Transient')
    else
        error('Heat transfer model not identified')
    end
    


    a_T(1,2) = -max(0, -rho_f(2,j)*v(2,j)*cp_f(2,j-1)/dx);    

    a_T(1,1) = rho(1,j)*cp(1,j-1)/dt ...
        + (rho_f(2,j)*v(2,j)*cp_f(2,j-1) - rho_f(1,j)*v(1,j)*cp_f(1,j-1))/dx...
        - (a_T(1,2) - max(rho_f(1,j)*v(1,j)*cp_f(1,j-1)/dx, 0));
    
    b_T(1) = Q(j) + (P(1,j)-P(1,j-1))/dt ...
            + v_n(1,j)*(P_f(2,j) - P_f(1,j))/dx ...
            + f(1)*rho(1,j)*abs(v_n(1,j))^3/(2*D) ...
            + rho(1,j-1)*cp(1,j-1)*T(1,j-1)/dt...
            + max(rho_f(1,j)*v(1,j)*cp_f(1,j-1)/dx, 0)*T_f(1,j);


    a_T(n_n,n_n-1) = -max(rho_f(n_f-1,j)*v(n_f-1,j)*cp_f(n_f-1,j-1)/dx, 0);

    a_T(n_n,n_n) = rho(n_n,j)*cp(n_n,j-1)/dt ...
        + (rho_f(n_f,j)*v(n_f,j)*cp_f(n_f,j-1) - rho_f(n_f-1,j)*v(n_f-1,j)*cp_f(n_f-1,j-1))/dx...
        - (a_T(n_n,n_n-1) - max(0, -rho_f(n_f,j)*v(n_f,j)*cp_f(n_f,j-1)/dx));
    
    b_T(n_n) = Q(j) + (P(n_n,j)-P(n_n,j-1))/dt ...
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
        
        b_T(i) = Q(j) + (P(i,j)-P(i,j-1))/dt ...
                + v_n(i,j)*(P_f(i+1,j) - P_f(i,j))/dx ...
                + f(i)*rho(i,j)*abs(v_n(i,j))^3/(2*D) ...
                + rho(i,j-1)*cp(i,j-1)*T(i,j-1)/dt;
    end
    
    T(:,j) = linsolve(a_T,b_T);
    
    % THERMODYNAMIC PROPERTIES

    for i = 1:n_n  % PROPERTIES FROM P AND RHO          
        % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),'Air');
        cp(i,j) = CP.PropsSI('C','P',P(i,j),'D',rho(i,j),'Air');
        h(i,j) = CP.PropsSI('H','P',P(i,j),'D',rho(i,j),'Air');
        s(i,j) = CP.PropsSI('S','P',P(i,j),'D',rho(i,j),'Air');
        u(i,j) = CP.PropsSI('U','P',P(i,j),'D',rho(i,j),'Air');
        cp_f(i,j) = CP.PropsSI('C','P',P_f(i,j),'D',rho_f(i,j),'Air');
        h_f(i,j) = CP.PropsSI('H','P',P_f(i,j),'D',rho_f(i,j),'Air');
        s_f(i,j) = CP.PropsSI('S','P',P_f(i,j),'D',rho_f(i,j),'Air');
    end
    
    cp_f(n_f,j) = CP.PropsSI('C','P',P_f(n_f,j),'D',rho_f(n_f,j),'Air');
    h_f(n_f,j) = CP.PropsSI('H','P',P_f(n_f,j),'D',rho_f(n_f,j),'Air');
    s_f(n_f,j) = CP.PropsSI('S','P',P_f(n_f,j),'D',rho_f(n_f,j),'Air');

    % Upwind scheme
    T_f(2:end-1,j) = (v(1:end-2,j) >= 0).*T(1:end-1,j) ...
                     + (v(1:end-2,j) <  0).*T(2:end,j);
    
    T_f(end,j) = T(end,j); % ASSUMING UPWIND SCHEME !!!!!!!!!!
                           % Implement other cases
    

    % CAES process
    if strcmp(Process,'Charging_L')
        if strcmp(L_bound,'Inlet') & P(ceil(n_n/2),j) >= P_max
        % if strcmp(L_bound,'Inlet') & P(1,j) >= P_max
            % Charging from L boundary
            % v(1,j+1:end) = 0;
            L_bound = 'Wall';
            t_shut_off = (j-1)*dt; 
        end
    elseif strcmp(Process,'Discharging_L')
        % if strcmp(L_bound,'M_const') & P(1,j) <= P_min
        if strcmp(L_bound,'M_const') & P(ceil(n_n/2),j) <= P_min
            % Discharging from L boundary
            L_bound = 'Wall';
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
        elseif strcmp(stage,'idle_charg') & t(j) >= Dt_charg
            L_bound = 'M_const';
            j_disch = j; % index of start of discharge
            m_L = -100;
            stage = 'Discharging';
        elseif strcmp(stage,'Discharging') & P(1,j) <= P_min
            L_bound = 'Wall';
            stage = 'idle_disch';
            
        end
        stage_hist = {stage_hist;stage};
    elseif strcmp(Process,'Cycle_R')
    elseif strcmp(Process,'NCycle_L')
    elseif strcmp(Process,'NCycle_R')
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
        disp('20% of simulation time.')
    elseif t(j-1) < 2*Dt/5 && t(j) >= 2*Dt/5
        disp('40% of simulation time.')
    elseif t(j-1) < 3*Dt/5 && t(j) >= 3*Dt/5
        disp('60% of simulation time.')
    elseif t(j-1) < 4*Dt/5 && t(j) >= 4*Dt/5
        disp('80% of simulation time.')
    end

end

elapsedTime = toc

dm = zeros(n_t,1);
dE = zeros(n_t,1);
% m_bal = zeros(n_t,1);
% E_bal = zeros(n_t,1);

% As a function of inlet velocity (CONSTANT INLET CONDITIONS)
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
else
    error('Process not identified.')
end

% m_bal = m(1) + cumsum([0;dm(1:end-1)]);
% E_bal = E(1) + cumsum([0;dE(1:end-1)]);
m_bal = m(1) + cumsum(dm);
E_bal = E(1) + cumsum(dE);

x = 0:dx:L;
x_f = [0:dx:L]';
x_n = [dx/2:dx:L]';

% Exergy
% X = P*(A_h*L).*(P_amb./P - 1 + log(P./P_amb))./(1e6*3600);          % Pipeline Exergy [MWh]
% X_min = P_0*(A_h*L).*(P_amb./P_0 - 1 + log(P_0./P_amb))/(1e6*3600); % Exergy when discharged [MWh] 
% X = P*(A_h*dx).*(P_amb./P - 1 + log(P./P_amb))./(1e6*3600);           % Pipeline Exergy [MWh]
% X_min = P_0*(A_h*dx).*(P_amb./P_0 - 1 + log(P_0./P_amb))/(1e6*3600);  % Exergy when discharged [MWh] 

X = sum(m_n.*( u - u_o + P_o*R*(T./P - T_o/P_o) - T_o*(s - s_o) ) )./(1e6*3600);           % Pipeline Exergy [MWh]

X_0 = m(1)*(u_0 - u_o + P_o*R*(T_0./P_0 - T_o/P_o) - T_o*(s_0 - s_o) )/(1e6*3600);  % Exergy at t = 0s [MWh] 


X_net = X - X_0;                                                 % Exergy between current state and discharged state (assuming whole pipeline at P_min)
DeltaX_flow = sum(dX)/(1e6*3600);
DeltaX_st = sum(X_net(:,end));
% error_X = (DeltaX_st - DeltaX_flow)/DeltaX_st

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

% 
% if strcmp(Process,'Charging_L') || strcmp(Process,'Charging_R')
% 
% elseif strcmp(Process,'Discharging_L') || strcmp(Process,'Discharging_R')
% 
% end

% Figures of mass and energy over time
% if legend('X','X_0 + dX')
%     figure('color',[1 1 1]);plot(m)
%     hold on; plot(mm)
%     legend('m','\Delta m')
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
else
    if Figures
        figure('color',[1 1 1]);plot(t,m)
        hold on; plot(t,m_bal)
        legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
        figure('color',[1 1 1]);plot(t,E./(1e6*3600))
        hold on; plot(t,E_bal(1:end)./(1e6*3600))
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
    yyaxis right
    plot(t,E2,t,E2_bal)
    legend('m','m_{bal}','E','E_{bal}')
    
    figure
    plot((E2 - E2_bal')./E2)
end

XX = sum(m_n(2:end,:).*( u(2:end,:) - u_o + P_o*R*(T(2:end,:)./P(2:end,:) - T_o/P_o) - T_o*(s(2:end,:) - s_o) ) )./(1e6*3600);

dXX = rho_f(2,:)'.*v(2,:)'*A_h.*(h_f(2,:)' - h_o - T_o*(s_f(2,:)' - s_o))*dt./(1e6*3600);

XX_bal = XX(1) + cumsum(dXX);

etaX_stor = XX(end)/XX_bal(end)


% filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'MPa_m_',num2str(m_in),'_',string(year(time)),'_',string(month(time)),'_',string(day(time)));
filename = strcat(simType,'_',Process,'_P_',num2str(P_max/1e6),'-',num2str(P_min/1e6),'MPa_m_',num2str(m_in));
save(filename)

L_bound
%% Previous tests

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


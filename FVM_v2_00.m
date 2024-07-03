clear 
% close all
clc

% profile on
tic
CP = py.importlib.import_module('CoolProp.CoolProp'); % Simplifies coolprop calls

%----------------------- PROBLEM PARAMETERS -------------------------%
% Pipeline parameters
% Parameters from Kiuchi (
L = 5000;
D = 0.5;
% L = 70000; % Reference case 70 km length, 0.9 m diam
% D = 0.9;
% L = 213333; % Length for eq. vol with Huntorf

% Cavern parameters
% L = 300; % Volume shape similar to Huntorf
% D = 24;
% L = 35; % Roughly same volume as 70 km, 0.9 D pipeline
% D = 40;

% Total simulation time
% Dt = 3*3600; 
% Dt = 3600;
Dt = 180;

eps = 0.04e-3; % Absolute roughness 0.04 mm

% Ambient conditions
P_amb = 101325;
T_amb = 273.15 + 25;

% System operational limits
P_max = 7e6;
P_min = 4.3e6;

% Type of simulation - Cavern or pipeline storage
simType = 'CAESPipe';
% simType = 'CAESCav';

% CAES process
Process = 'Charging_L';
% Process = 'Discharging_L'
% Process = 'Charging_R';
% Process = 'Discharging_R';
if strcmp(Process,'Charging_L')
    % Initial conditions
    % P_0 = 101325;
    P_0 = 4.3e6; % Huntorf
    T_0 = 273.15 + 25;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Inlet'; % OR M_CONST ??
    R_bound = 'Wall';
elseif strcmp(Process,'Discharging_L')
    % Initial conditions
    P_0 = 7e6;
    T_0 = 273.15 + 25;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'M_const';
    R_bound = 'Wall';
elseif strcmp(Process,'Charging_R')
    % Initial conditions
    % P_0 = 101325;
    P_0 = 4.3e6; % Huntorf
    T_0 = 273.15 + 25;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Wall';
    R_bound = 'Inlet';
    
elseif strcmp(Process,'Discharging_R')
    % Initial conditions
    P_0 = 7e6;
    T_0 = 273.15 + 25;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Wall';
    R_bound = 'M_const';
end

% A - Left side boundary condition
% L_bound = 'Outlet';
% L_bound = 'Wall';
% L_bound = 'Inlet';
% L_bound = 'P_const';
% L_bound = 'M_const';
if strcmp(L_bound,'Inlet')
    % Q_st_in = 3e5; % standard cubic meters per hour
    % Q_a = Q_st_in/3600;
    % m_in = rho_a*Q_a;
    % m_in = 100;
    m_L = 108; % Huntorf
    % v_A = ?
    P_L = P_0;
    T_L = 273.15 + 25;

elseif strcmp(L_bound,'Wall')
    m_L = 0;
    v_L = 0;
elseif strcmp(L_bound,'Outlet')
    % m_A = m_R;
elseif strcmp(L_bound,'P_const')
    P_L = 4.3e6; % 4 MPa
elseif strcmp(L_bound,'M_const')
    % m_A = -417; % Huntorf, sign indicates flow direction
    m_L = -100;
end   

% B - Right side boundary condition
% R_bound = 'Outlet';
% R_bound = 'Inlet';
% R_bound = 'P_const';
% R_bound = 'Wall';
% R_bound = 'M_const';
if strcmp(R_bound,'Outlet')
    % m_out = m_A;
elseif(strcmp(R_bound,'Wall'))
    m_R = 0;
    v_R = 0;
elseif(strcmp(R_bound,'Inlet'))
    m_R = -108; % Huntorf  % Remember the sign to indicate the direction of flow (+ right , - left)
    P_R = 7e6;
    T_R = 273.15 + 60;
    % v_B = m_B/(rho_B*A_h);
elseif strcmp(R_bound,'P_const')
    P_R = 4e6; 
    % P_B = 4.13e6; % 4.13 MPa = Huntorf https://www.sciencedirect.com/science/article/pii/S0196890420302004#s0010
elseif strcmp(R_bound,'M_const')
    % m_B = 417; % Huntorf    
    m_R = 100;
end                     

g = 9.81;
theta = 0;

%--------------------- SIMULATION PARAMETERS ------------------------%

if strcmp(simType,'CAESPipe')
    dx = L/(40-1); % CAESPipe
    dt = 1;
    tol = 1e-6; % CAESPipe Charging
    % tol = 1e-3; % CAESPipe discharging
elseif strcmp(simType,'CAESCav')
    dx = L/(5-1); % CAESCav
    dt = 0.01;
    tol = 1e-7; % CAEScav
else
    warning('Cant identify simulation type.')
end

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
% h_f = zeros(n_f,n_t);
% s_f = zeros(n_f,n_t);
% u_f = zeros(n_f,n_t);
% cp_f = zeros(n_f,n_t);

% Node initialization
v_n = zeros(n_n,n_t);
P = zeros(n_n,n_t);
T = zeros(n_n,n_t);
rho = zeros(n_n,n_t);
h = zeros(n_n,n_t);
s = zeros(n_n,n_t);
u = zeros(n_n,n_t);
cp = zeros(n_n,n_t);
T_sol = zeros(n_n,n_t); % temporary var to store T as solved using energy 
% equation and compare to using only mass and momentum eqs

cp_f = zeros(1,n_t);
h_f = zeros(1,n_t);
s_f = zeros(1,n_t);

% Sanity check
m = zeros(n_t,1); % Total mass in pipeline
E = zeros(n_t,1); % Total energy in pipeline

m_n = zeros(n_n,n_t); % Mass per node
E_n = zeros(n_n,n_t); % Energy per node

%------------------ SETTING INITIAL CONDITIONS ---------------------%
% Basic calculations
epsD = eps/D;
A_h = pi*D^2/4;

rho_amb = CP.PropsSI('D','P',P_amb,'T',T_amb,'Air');
% Dead state parameters
P_o = P_amb;
T_o = T_amb;
h_o = CP.PropsSI('H','P',P_o,'T',T_o,'Air');
s_o = CP.PropsSI('S','P',P_o,'T',T_o,'Air');

% Initial conditions at nodes (i -> x, j -> t)
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

P = 1:10;
v = -5:5;

% Properties in nodes to faces
P_f(2:end-1,1) = (v(2:end-1,1) >= 0).*P(1:end-1,1) ...
    +            (v(2:end-1,1) <  0).*P(2:end,1);
T_f(2:end-1,1) = (v(2:end-1,1) >= 0).*T(1:end-1,1) ...
    +            (v(2:end-1,1) <  0).*T(2:end,1);
rho_f(2:end-1,1) = (v(2:end-1,1) >= 0).*rho(1:end-1,1) ...
    +            (v(2:end-1,1) <  0).*rho(2:end,1);
% 
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

    % v_L = m_L/(rho_L*A_h);

    % P(1,:) = P_L;
    % T(1,:) = T_A;
    % rho(1,:) = rho_A;

    P_f(1,1) = P_L;
    T_f(1,:) = T_L;
    rho_f(1,1) = rho_L;
    
    % v_L = m_L/(rho_f(1,1)*A_h);
    v(1,1) = m_L/(rho_f(1,1)*A_h);

elseif strcmp(L_bound,'Wall')
    v(1,:) = 0;

    P_f(1,1) = P(1,1);
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
    P(1,1) = P(2,1);
    T(1,1) = T(2,1);
    rho(1,1) = rho(2,1);

    v_n(1,1) = m_L/(rho(1,1)*A_h);

    % Upwind scheme
    v(1,1) = v_n(1,1);
    
    P_f(1,1) = P(1,1);
    T_f(1,1) = T(1,1);
    rho_f(1,1) = rho(1,1);
    
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
    v(end,1) = m_R/(rho_f(end,1)*A_h);

elseif(strcmp(R_bound,'Wall'))
    v(end,:) = 0;

    P_f(end,1) = P(end,1);
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

m(1) = sum(rho(:,1)*A_h*dx);
E(1) = sum(rho(:,1)*A_h*dx.*u(:,1));

m_n(:,1) = rho(:,1)*A_h*dx;
E_n(:,1) = rho(:,1)*A_h*dx.*u(:,1);

count = zeros(n_t,1);

% Friction factor based on Nikuradse
f_guess = (2*log10(1/epsD)+1.14)^(-2);


error_hist = [];
bound_hist = [string(L_bound) string(R_bound)];

for j=2:n_t
    % Initial guess for next time step is the same props as the previous t step
    


    % Momentum and mass balance loop
    


end
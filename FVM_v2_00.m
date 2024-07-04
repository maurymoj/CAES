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
    % P_R = 7e6;
    P_R = P_0; % Sliding pressure at pipeline inlet
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

% SHOULD THESE BE MATRICES OR VECTORS ???
cp_f = zeros(1,n_t);
h_f = zeros(1,n_t);
s_f = zeros(1,n_t);

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
    v(1,1) = m_L/(rho_f(1,1)*A_h);

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
    P(:,j) = P(:,j-1);
    rho(:,j) = rho(:,j-1);
    T(:,j) = T(:,j-1);



    v(:,j) = v(:,j-1);
    % v(2:end-1,j) = v(2:end-1,j-1);


    P_corr = zeros(n_n,1);
    rho_corr = zeros(n_n,1);
    v_corr = zeros(n_f,1);
    
    v_star = v(:,j);
    count(j) = 0;
    error_P = 10;

    % Momentum and mass balance loop
    while count(j) < 100 && max(abs(error_P)) > tol
        % Under-relaxed corrections
        P(:,j) = P(:,j) + alpha_P*P_corr;
        rho(:,j) = alpha_rho*(rho(:,j) + rho_corr) + (1-alpha_rho)*rho(:,j);

        v(:,j) = alpha_v*(v_star + v_corr) + (1-alpha_v)*v_star;
        % v(1:end-1,j) = alpha_v*(v_star(1:end-1) + v_corr(1:end-1)) + (1-alpha_v)*v_star(1:end-1);

        % Boundary conditions
        if strcmp(L_bound,'Inlet')
            % P(1,:) = P_L;
            % T(1,:) = T_L;
            % rho(1,:) = rho_L;

            P_f(1,j) = P(1,j); % Assumption of constant pressure between 
                               % face and first node
            T_f(1,j) = T_L; 
            % T_f(1,j) = T(1,j); 
            rho_f(1,j) = CP.PropsSI('D','P',P_f(1,j),'T',T_f(1,j),'Air');
            v(1,j) = m_L/(rho_f(1,j)*A_h); % Sliding pressure test
            
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
            % Assumptions:
            % - Zero gradient for all properties except u-velocity
            % - Constant cross-sectional area
            P(1,j) = P(2,j);
            T(1,j) = T(2,j);
            rho(1,j) = rho(2,j);
        
            v_n(1,j) = m_L/(rho(1,j)*A_h);
        
            % Upwind scheme
            v(1,j) = v_n(1,j);
            
            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);
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
            
            % v_n(end,j) = v(end-1,j);
            v_n(end,j) = (v(end-1,j) >= 0).*v(end-1,j) ...
            +      (v(end-1,j) <  0).*v(end,j); 
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
            nu = CP.PropsSI('viscosity','P',P_f(i,j),'D',rho_f(i,j),'Air');

            drho_dP_f(i) = 1/(u_sonic_f(i)^2);

            u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),'Air');
            drho_dP_n(i) = 1/(u_sonic_n(i)^2);
        end

        u_sonic_f(end) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'D',rho_f(end,j),'Air');
        drho_dP_f(end) = 1/(u_sonic_f(end)^2);



        %---------- v* calculation - Momentum control volume ---------------------%
        
        % Friction factor based on Nikuradse - IMPLEMENT COLEBROOK EQUATION
        f = (2*log10(1/epsD)+1.14)^(-2)*ones(n_n,1);
        
        % Re = rho_f(:,j).*v(:,j)*D./nu(:);
        % f = zeros(n,1);
        % 
        % for i=1:n_n
        %     f_old = f_guess;
        %     df = 10;
        %     count_f = 0;
        %     while (df > 0.0001 & count_f < 10) 
        %         % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
        %         f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
        %         df = abs((f_new - f_old)/f_old);
        %         f_old = f_new;
        %         count_f = count_f + 1; 
        %     end
        %     f(i) = f_old;
        % end
        


        % Matrix/vector initialization (only need to solve for the inner faces -
        % boundary faces are dealt with the boundary conditions)
        a_M = zeros(n_f);
        b_M = zeros(n_f,1);
        d_M = zeros(n_f,1);
        B_M = zeros(n_f,1);
        
        a_M(1,1) = 1; % Value for v at the first momentum volume is 
                    % known from boundary conditions
        B_M(1) = v(1,j);
        
        a_M(end,end) = 1;
        B_M(end) = v(end,j);
        

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
        % WHAT TO DO REGARDING TO THE a INDICES ??
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
        
        % v_star = B\a
        v_star = linsolve(a_M,B_M);
        % v_star = [v(1,j);v_star;v(end,j)];
        % if strcmp(L_bound,'Outlet')
        %     v_star(1) = v_star(2);
        % elseif strcmp(R_bound,'Outlet')
        %     v_star(end) = v_star(end-1);
        % end

    end


end
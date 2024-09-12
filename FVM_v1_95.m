clear 
% close all
% clc
% profile on
tic
CP = py.importlib.import_module('CoolProp.CoolProp');

%----------------------- PROBLEM PARAMETERS -------------------------%
% Pipeline parameters
% Kiuchi
L = 5000;
D = 0.5;
% L = 70000;
% D = 0.9;
% L = 213333;
% Cavern parameters
% L = 300;
% D = 24;
% L = 35; % Roughly same volume as 70 km, 0.9 D pipeline
% D = 40;

% Dt = 3*3600; % Total simulation time
% Dt = 3600;
Dt = 90;

eps = 0.04e-3; % Absolute roughness 0.04 mm

% Ambient conditions
P_amb = 101325;
T_amb = 273.15 + 25;

% System operational limits
P_max = 7e6;
P_min = 4e6;

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
    P_0 = 4e6; % Huntorf
    T_0 = 273.15 + 25;
    v_0 = 0;
    % v_0 = v_in;

    L_bound = 'Inlet';
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
    T_L = 273.15 + 60;

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
N_f = n_n + 1;      % number of faces
n_t = Dt/dt+1;  % n of time steps

% Face initialization
v = zeros(N_f,n_t);
P_f = zeros(N_f,n_t);
T_f = zeros(N_f,n_t);
rho_f = zeros(N_f,n_t);

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

% Properties in nodes to faces
P_f(2:end-1,1) = (v(1:end-2,1) >= 0).*P(1:end-1,1) ...
    +            (v(1:end-2,1) <  0).*P(2:end,1);
T_f(2:end-1,1) = (v(1:end-2,1) >= 0).*T(1:end-1,1) ...
    +            (v(1:end-2,1) <  0).*T(2:end,1);
rho_f(2:end-1,1) = (v(1:end-2,1) >= 0).*rho(1:end-1,1) ...
    +            (v(1:end-2,1) <  0).*rho(2:end,1);

P_f(end,1) = P(end,1); % ASSUMING v >= 0 for t=0!!!!
T_f(end,1) = T(end,1); % ASSUMING v >= 0 for t=0!!!!
rho_f(end,1) = rho(end,1); % ASSUMING v >= 0 for t=0!!!!


% Boundary conditions (Inlet, Wall, and constant pressure)

% L boundary conditions
if strcmp(L_bound,'Inlet')
    % At node
    rho_L = CP.PropsSI('D','P',P_L,'T',T_L,'Air');
    cp_L = CP.PropsSI('C','P',P_L,'T',T_L,'Air');
    h_L = CP.PropsSI('H','P',P_L,'T',T_L,'Air');
    s_L = CP.PropsSI('S','P',P_L,'T',T_L,'Air');

    % v_L = m_L/(rho_L*A_h);

    % P(1,:) = P_L;
    % T(1,:) = T_A;
    % rho(1,:) = rho_A;

    P_f(1,:) = P_L;
    T_f(1,:) = T_L;
    rho_f(1,:) = rho_L;
    % At face
    % v(1,2:end) = v_L;
    
    % Sliding pressure inlet test

    v_L = m_L/(rho(1,1)*A_h);
    v_n(1,1) = 0;
    % P(1,:) = P_0;
    % T(1,:) = T_A;
    % rho(1,:) = rho_A;
    % 
    % P_f(1,:) = P_0;
    % T_f(1,:) = T_0;
    % rho_f(1,:) = rho(1,1);
    % % At face
    % v(1,1) = v_A;

elseif strcmp(L_bound,'Wall')
    v(1,:) = 0;

    % P_f(1,1) = P(1,1);
    % T_f(1,1) = T(1,1);
    % rho_f(1,1) = rho(1,1);

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

    v_R = m_R/(rho_R*A_h);

    P(end,:) = P_R;
    T(end,:) = T_R;
    rho(end,:) = rho_R;

    P_f(end,:) = P_R;
    T_f(end,:) = T_R;
    rho_f(end,:) = rho_R;
    % At face
    v(end,:) = v_R;

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

f_guess = (2*log10(1/epsD)+1.14)^(-2); % Friction factor based on Nikuradse



error_hist = [];
bound_hist = [string(L_bound) string(R_bound)];

for j=2:n_t
    % Initial guess for next time step is the same props as the previous t step
    P(:,j) = P(:,j-1);
    rho(:,j) = rho(:,j-1);
    T(:,j) = T(:,j-1);
    v(2:end-1,j) = v(2:end-1,j-1);
    
    P_corr = zeros(n_n,1);
    rho_corr = zeros(n_n,1);
    v_corr = zeros(N_f,1);
    
    v_star = v(:,j);
    count(j) = 0;
    error_P = 10;

    while count(j) < 100 && max(abs(error_P)) > tol
        
        % Under-relaxed corrections
        P(:,j) = P(:,j) + alpha_P*P_corr;
        rho(:,j) = alpha_rho*(rho(:,j) + rho_corr) + (1-alpha_rho)*rho(:,j);

        v(1:end-1,j) = alpha_v*(v_star(1:end-1) + v_corr(1:end-1)) + (1-alpha_v)*v_star(1:end-1);
        
        % L boundary
        if strcmp(L_bound,'Inlet')
            % P(1,:) = P_L;
            % T(1,:) = T_L;
            % rho(1,:) = rho_L;

            P_f(1,j) = P(1,j);
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
        elseif strcmp(L_bound,'P_const')
            % P(1,:) = P_L;
            T(1,j) = T(2,j);
            rho(1,j) = rho(2,j);

            v(1,j) = v(2,j);
        
            P_f(1,j) = P(1,j);
            T_f(1,j) = T(1,j);
            rho_f(1,j) = rho(1,j);
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
        elseif(strcmp(R_bound,'Wall'))
            % Set by upwind scheme
            v(end,j) = 0;
        elseif strcmp(R_bound,'P_const')
            % P(end,:) = P_R;
            T(end,j) = T(end-1,j);
            rho(end,j) = rho(end-1,j);
        
            v(end,j) = v(end-1,j);
        
            P_f(end,j) = P(end,j);
            T_f(end,j) = T(end,j);
            rho_f(end,j) = rho(end,j);   
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

        % Initial guess (properties at t+dt = properties at t) - without
        % under-relaxation
        % P(:,j) = P(:,j) + P_corr;
        % rho(:,j) = rho(:,j) + rho_corr;
        % T(:,j) = T(:,j);
        % v(:,j) = v_star(:) + v_corr;
        
        % Upwind scheme
        % v_n(:,j) = (v(1:end-1,j) >= 0).*v((1:end-1),j) ...
        %     +      (v(1:end-1,j) <  0).*v((2:end),j);
        v_n(2:end,j) = (v(2:end-1,j) >= 0).*v((2:end-1),j) ...
            +      (v(2:end-1,j) <  0).*v((3:end),j);
        
        % Properties in nodes to faces
        P_f(2:end-1,j) = (v(1:end-2,j) >= 0).*P(1:end-1,j) ...
            +            (v(1:end-2,j) <  0).*P(2:end,j);
        rho_f(2:end-1,j) = (v(1:end-2,j) >= 0).*rho(1:end-1,j) ...
            +            (v(1:end-2,j) <  0).*rho(2:end,j);

        P_f(end,j) = P(end,j); % ASSUMING v >= 0 for t>0!!!!
        rho_f(end,j) = rho(end,j); % ASSUMING v >= 0 for t>0!!!!

        

        % Properties
        u_sonic = zeros(N_f,1);
        drho_dP = zeros(N_f,1);
        nu = zeros(N_f,1);

        u_sonic_n = zeros(n_n,1);
        drho_dP_n = zeros(n_n,1);

        for i = 1:n_n  % PROPERTIES FROM P AND RHO
            % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),'Air');

            u_sonic(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),'Air');
            nu = CP.PropsSI('viscosity','P',P_f(i,j),'D',rho_f(i,j),'Air');

            drho_dP(i) = 1/(u_sonic(i)^2);

            u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),'Air');
            drho_dP_n(i) = 1/(u_sonic_n(i)^2);
        end

        u_sonic(end) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'D',rho_f(end,j),'Air');
        drho_dP(end) = 1/(u_sonic(end)^2);

        % T_f(2:end-1,j) = (v(1:end-2,j) >= 0).*T(1:end-1,j) ...
        %     +            (v(1:end-2,j) <  0).*T(2:end,j);
        % T_f(end,j) = T(end,j); % ASSUMING v >= 0 for t>0!!!!
        % 
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
        
        % Re = rho*
        % f_old = f_guess;
        % df = 10;
        % count_f=0;
        % while (df > 0.0001 & count_f < 10) 
        %     % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
        %     f_new = (-2*log10(epsD/3.7 + 2.825/(Re*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
        %     df = abs((f_new - f_old)/f_old);
        %     f_old = f_new;
        %     count_f = count_f + 1; 
        % end
        % f = f_old;


        % Matrix/vector initialization (only need to solve for the inner faces -
        % boundary faces are dealt with the boundary conditions)
        a = zeros(N_f-2);
        b = zeros(N_f-2,1);
        d = zeros(N_f-2,1);
        B = zeros(N_f-2,1);
        
        % For momentum equations the first line refers to face B and the last line to
        % face E, as we know velocities at faces A and F from the boundary
        % conditions
        a(1,2) = -max( 0, -rho(2,j)*v_n(2,j)/dx);    % a_C
        a(1,1) = rho_f(2,j)/dt + f(1)*rho_f(2,j)*abs(v(2,j))/(2*D) ...
            + (rho(2,j)*v_n(2,j) - rho(1,j)*v_n(1,j))/dx ...
            - (a(1,2) - max(rho(1,j)*v_n(1,j)/dx,  0) );% a_B
                
        d(1) = -1/dx;
        
        
        
        
        
        
        
        % b(1) = (rho_f(2,j-1)*v(2,j-1)/dt)...
        %     - rho_f(2,j)*g*sind(theta)...
        %     + max(rho(1,j)*v_n(1,j)/dx,  0)*v_n(1,j);
        b(1) = (rho_f(2,j-1)*v(2,j-1)/dt)...
            - rho_f(2,j)*g*sind(theta)...
            + max(rho(1,j)*v_n(1,j)/dx,  0)*v_n(1,j);





        
        B(1) = d(1)*(P(2,j)-P(1,j)) ...
            + b(1);           
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
        
        a(end,end-1) = -max(rho(n_n-1,j)*v_n(n_n-1,j)/dx, 0); % a_N-2 (index N-2, N-3)
        % WHAT TO DO REGARDING TO THE a INDICES ??
        a(end,end) = ...                                  % a_N-1 (index N-2, N-2)
            rho_f(N_f-1,j)/dt + f(N_f-1)*rho_f(N_f-1,j)*abs(v(N_f-1,j))/(2*D)...
            +(rho(n_n,j)*v_n(n_n,j) - rho(n_n-1,j)*v_n(n_n-1,j))/dx...
            - (a(end,end-1) - max(-rho(n_n,j)*v_n(n_n,j)/dx,  0) );
        
        d(end) = -1/dx;
        b(end) = (rho_f(N_f-1,j-1)*v(N_f-1,j-1)/dt)-rho_f(N_f-1,j)*g*sind(theta);
        
        B(end) = d(end)*(P(n_n,j)-P(n_n-1,j)) + b(end);
        
        for i=2:N_f-3
        % The indices of the node and face properties have to be added by 1, since
        % the 1st line was ommited (as we have it solved from the boundary
        % condition 
            a(i,i-1) = -max(rho(i,j)*v_n(i,j)/dx,                     0);
            a(i,i+1) = -max(                       0, -rho(i+1,j)*v_n(i+1,j)/dx);
            a(i,i)   = rho_f(i+1,j)/dt + f(i+1)*rho_f(i+1,j)*abs(v(i+1,j))/(2*D) ...
                - (a(i,i-1) + a(i,i+1)) ...
                + (rho(i+1,j)*v_n(i+1,j)/dx - rho(i,j)*v_n(i,j)/dx);
            d(i) = -1/dx;
            b(i) = (rho_f(i+1,j-1)*v(i+1,j-1)/dt)-rho_f(i+1,j)*g*sind(theta);
            
            B(i) = d(i)*(P(i+1,j)-P(i,j)) + b(i);
        end
        
        % v_star = B\a
        v_star = linsolve(a,B);
        v_star = [v(1,j);v_star;v(end,j)];
        % if strcmp(L_bound,'Outlet')
        %     v_star(1) = v_star(2);
        % elseif strcmp(R_bound,'Outlet')
        %     v_star(end) = v_star(end-1);
        % end
        % v(:,j) = [v(1,j);v_star;v(end,j)];

        %--------------- Pressure correction P' calculations ---------------------%
        % The matrices/vectors in P' calculations have no shift, so index 1 means
        % control volume 1 and so on
        A = zeros(n_n);
        BB = zeros(n_n,1);
        
        % left node
        % Coefficients of a and d are displaced by 1 (1 = B, 2 = C and so on) 
        % A(1,2) = - max(0, -drho_dP(2)*v_star(2)/dx); % A_2
        % A(1,1) = drho_dP_n(1)/dt - (rho_f(2,j)*d(1)/a(1,1))/dx ...
        %        + drho_dP(2)*v_star(2)/dx - A(1,2); % A_1
        A(1,2) = (rho_f(2,j)*d(1)/a(1,1))/dx - max(0, -drho_dP(2)*v_star(2)/dx); % A_2
        A(1,1) = drho_dP_n(1)/dt + drho_dP(2)*v_star(2)/dx ...
               - A(1,2); % A_1
        
        BB(1) = (rho(1,j-1)-rho(1,j))/dt + (rho_f(1,j)*v_star(1) - rho_f(2,j)*v_star(2))/dx;
        if strcmp(L_bound,'P_const')
            A(1,1) = 1;
            A(1,2) = 0;

            B(1) = 0;
        elseif strcmp(L_bound,'Outlet')
            %???
            A(1,2) = -drho_dP(2)*v_star(2)/dx + rho_f(2,j)*(d(1)/a(1,1))/dx; % A_2
            A(1,1) = drho_dP_n(1)/dt ...
                     - (rho_f(2,j)*d(1)/a(1,1))/dx ...
                     - drho_dP(1)*v_star(1)/dx ; % A_1
        end
        
        % right node
        A(n_n,n_n-1) = (rho_f(N_f-1,j)*d(end)/a(end,end))/dx...
                       - max(drho_dP(N_f-1)*v_star(N_f-1)/dx, 0);
        A(n_n,n_n) = drho_dP_n(n_n)/dt ...
                     - drho_dP(N_f-1)*v_star(N_f-1)/dx... 
                     - A(n_n,n_n-1); % v_star(end) = 0
        
        BB(n_n) = (rho(n_n,j-1)-rho(n_n,j))/dt + (rho_f(N_f-1,j)*v_star(N_f-1) - rho_f(N_f,j)*v_star(N_f))/dx;
        if strcmp(R_bound,'P_const')
            A(n_n,n_n-1) = 0;
            A(n_n,n_n) = 1;

            B(n_n) = 0;
        elseif strcmp(R_bound,'Outlet')
            %???
            A(n_n,n_n-1) = (rho_f(N_f-1,j)*d(end)/a(end,end))/dx...
                           - drho_dP(N_f-1)*v_star(N_f-1)/dx;
            A(n_n,n_n) = drho_dP_n(n_n)/dt ...
                         - (rho_f(N_f-1,j)*d(end)/a(end,end))/dx ...
                         + drho_dP(N_f)*v_star(N_f)/dx;
        end

        for i=2:n_n-1
        % Subtracted 1 from a and d as there is no row for the first face -
        % boundary conditions for matrix/vector in momentum eq.
            A(i,i-1) = (rho_f(i,j)*d(i-1)/a(i-1,i-1))/dx ...
                       - max(drho_dP(i)*v_star(i)/dx,                     0);
            A(i,i+1) = (rho_f(i+1,j)*d(i)/a(i,i))/dx ...
                       - max(              0, -drho_dP(i+1)*v_star(i+1)/dx);
            A(i,i)   = drho_dP_n(i)/dt + (drho_dP(i+1)*v_star(i+1)/dx ...
                       - drho_dP(i)*v_star(i)/dx) - ( A(i,i-1) + A(i,i+1) );
            
            BB(i) = (rho(i,j-1)-rho(i,j))/dt ...
                    + (rho_f(i,j)*v_star(i) - rho_f(i+1,j)*v_star(i+1))/dx;
        end

        P_corr = linsolve(A,BB);
        
        if strcmp(L_bound,'M_const') %|| strcmp(L_bound,'Inlet')
            P_corr(1) = 0;
        elseif strcmp(R_bound,'M_const')
            P_corr(end) = 0;
        end

        rho_corr = drho_dP_n.*P_corr;

        a_i = diag(a);


        v_corr = zeros(N_f,1);
        % v_corr(1)         = d(1)./a_i(1).*P_corr(1);
        v_corr(2:end-1)   = d./a_i.*(P_corr(2:end)-P_corr(1:end-1));
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Inlet boundary condition
        v_corr(1)         = 0; 
        % Outlet boundary condition
        v_corr(end)       = 0;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        % v_corr(end)       =  % Velocity corrected based on mass balance (m_in =
        % m_out)
        error_P = (P_corr./P(:,j));
        error_hist = [error_hist error_P];

        % error_rho = (rho_corr./rho(:,j));
        % error_v = (v_corr./v(:,j));
        count(j) = count(j)+1;  

    end
    
    % ENERGY BALANCE
    % ASSUMING v>0
    a_T = zeros(n_n,n_n);
    b_T = zeros(n_n,1);
    Q(j) = 0; % ADIABATIC ASSUMPTION

    a_T(1,1) = rho(1,j)*cp(1,j-1)/dt + rho_f(2,j)*v(2,j)*cp(2,j-1)/dx;
    b_T(1) = Q(j) + (P(1,j)-P(1,j-1))/dt ...
            + v_n(1,j)*(P_f(2,j) - P_f(1,j))/dx ...
            + f(1)*rho(1,j)*abs(v_n(1,j))^3/(2*D) ...
            + rho(1,j-1)*cp(1,j-1)*T(1,j-1)/dt...
            + rho_f(1,j)*v(1,j)*cp(1,j-1)*T_f(1,j)/dx;
    for i=2:n_n
        a_T(i,i) = rho(i,j)*cp(i,j-1)/dt + rho_f(i+1,j)*v(i+1,j)*cp(i,j-1)/dx;
        a_T(i,i-1) = -rho_f(i,j)*v(i,j)*cp(i,j-1)/dx;
        b_T(i) = Q(j) + (P(i,j)-P(i,j-1))/dt ...
                + v_n(i,j)*(P_f(i+1,j) - P_f(i,j))/dx ...
                + f(i)*rho(i,j)*abs(v_n(i,j))^3/(2*D) ...
                + rho(i,j-1)*cp(i,j-1)*T(i,j-1)/dt;

        % a_T(i,i+1) = 
    end
    % T_sol(:,j) = linsolve(a_T,b_T);
    T(:,j) = linsolve(a_T,b_T);

    % THERMODYNAMIC PROPERTIES
    % cp(i,j) = CP.PropsSI('C','P',P(floor(end/2),j),'D',rho(floor(end/2),j),'Air');
    cp_f(j) = CP.PropsSI('C','P',P_f(1,j),'D',rho_f(1,j),'Air');
    h_f(j) = CP.PropsSI('H','P',P_f(1,j),'D',rho_f(1,j),'Air');
    s_f(j) = CP.PropsSI('S','P',P_f(1,j),'D',rho_f(1,j),'Air');

    for i = 1:n_n  % PROPERTIES FROM P AND RHO          
        % T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),'Air');
        cp(i,j) = CP.PropsSI('C','P',P(i,j),'D',rho(i,j),'Air');
        h(i,j) = CP.PropsSI('H','P',P(i,j),'D',rho(i,j),'Air');
        s(i,j) = CP.PropsSI('S','P',P(i,j),'D',rho(i,j),'Air');
        u(i,j) = CP.PropsSI('U','P',P(i,j),'D',rho(i,j),'Air');
    end

    T_f(2:end-1,j) = (v(1:end-2,j) >= 0).*T(1:end-1,j) ...
                     + (v(1:end-2,j) <  0).*T(2:end,j);
    T_f(end,j) = T(end,j);

    if strcmp(L_bound,'Inlet') & P(ceil(n_n/2),j) >= P_max
        % v(1,j+1:end) = 0;
        L_bound = 'Wall';
        t_shut_off = (j-1)*dt;
    elseif strcmp(L_bound,'M_const') & P(1,j) <= P_min
        L_bound = 'Wall';
    elseif strcmp(R_bound,'Inlet') & P(end-1,j) >= P_max
        % v(1,j+1:end) = 0;
        R_bound = 'Wall';
        t_shut_off = (j-1)*dt;
    elseif strcmp(R_bound,'M_const') & P(end,j) <= P_min
        R_bound = 'Wall';
    end

    bound_hist = [bound_hist; string(L_bound) string(R_bound)];
    
    m(j) = sum(rho(:,j)*A_h*dx);
    % E(j) = sum(rho(:,j)*A_h*dx.*cp(:,j).*T(:,j));
    E(j) = sum(rho(:,j)*A_h*dx.*u(:,j));
    
    m_n(:,j) = rho(:,j)*A_h*dx;
    E_n(:,j) = rho(:,j)*A_h*dx.*u(:,j);

    % Simulation iteration update
    % [j count(j)]
    % if rem((j-1)*dt,10) == 0
    %     disp(strcat('t= ',num2str((j-1)*dt),'s, n iterations: ',num2str(count(j))))
    % end
end

toc

% Sanity checks
% dm = rho_in*V_in*A_h*dt;
% dE = rho_in*V_in*A_h*cp_in*T_in*dt;
% mm = m(1) + [0; cumsum(dm(1:end-1))];
% EE = E(1) + [0; cumsum(dE(1:end-1))];
% mm = m(1) + dm.*(v~=0).*(0:n_t-1)';
% EE = E(1) + dE.*(0:n_t-1)';
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
    dE = rho_f(1,:)'.*v(1,:)'.*A_h.*h_f(:)*dt;
    dX = rho_f(1,:)'.*v(1,:)'*A_h.*(h_f(:) - h_o - T_o*(s_f(:) - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible
elseif strcmp(Process,'Discharging_L')
    dm = rho_f(1,:)'.*v(1,:)'*A_h*dt;
    dE = rho_f(1,:)'.*v(1,:)'*A_h.*cp(1,:)'.*T_f(1,:)'*dt;
    dX = rho_f(1,:)'.*v(1,:)'*A_h.*(h(1,:)' - h_o - T_o*(s(1,:)' - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible
elseif strcmp(Process,'Charging_R')
    dm = -rho_R*v(end,:)'*A_h*dt;
    dE = -rho_R*v(end,:)'*A_h*cp_R*T_R*dt;
    dX = -rho_R*v(end,:)'*A_h*(h_R - h_o - T_o*(s_R - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible
elseif strcmp(Process,'Discharging_R')
    dm = -rho_f(end,:)'.*v(end,:)'*A_h*dt;
    dE = -rho_f(end,:)'.*v(end,:)'*A_h.*cp(end,:)'.*T_f(end,:)'*dt;
    dX = -rho_f(end,:)'.*v(end,:)'*A_h.*(h(:,end)' - h_o - T_o*(s(:,end)' - s_o))*dt; % Flow exergy - kinetic and potential term contributions assumed negligible
end

m_bal = m(1) + cumsum(dm);
E_bal = E(1) + cumsum(dE);

x = 0:dx:L;
x_f = [0:dx:L]';
x_n = [dx/2:dx:L]';

% Exergy
% X = P*(A_h*L).*(P_amb./P - 1 + log(P./P_amb))./(1e6*3600);          % Pipeline Exergy [MWh]
% X_min = P_0*(A_h*L).*(P_amb./P_0 - 1 + log(P_0./P_amb))/(1e6*3600); % Exergy when discharged [MWh] 
X = P*(A_h*dx).*(P_amb./P - 1 + log(P./P_amb))./(1e6*3600);           % Pipeline Exergy [MWh]
X_min = P_0*(A_h*dx).*(P_amb./P_0 - 1 + log(P_0./P_amb))/(1e6*3600);  % Exergy when discharged [MWh] 

% IS THE e - eo TERM MISSING ??

X_net = X - X_min;                                                   % Exergy between current state and discharged state (assuming whole pipeline at P_min)
X_in = sum(dX)/(1e6*3600)
X_st = sum(X_net(:,end))

% Figures of mass and energy over time
% figure('color',[1 1 1]);plot(m)
% hold on; plot(mm)
% legend('m','\Delta m')

figure('color',[1 1 1]);plot(t,m)
hold on; plot(t,m_bal(1:end))
legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
figure('color',[1 1 1]);plot(t,E)
hold on; plot(t,E_bal(1:end))
legend('E','$E_o + \dot{m} \Delta E$','Interpreter','latex')
%%
% differences between total mass/energy in and change in C.V. mass/energy
% figure('color',[1 1 1]);plot(t,mm' - m)
% title('Difference in mass')
% figure('color',[1 1 1]);plot(t,EE' - E)
% title('Difference in energy')
% % Relative differences between total mass/energy in and change in C.V. mass/energy
figure('color',[1 1 1]);plot(t,(m_bal - m)./m)
title('Difference in mass')
figure('color',[1 1 1]);plot(t,(E_bal - E)./E)
title('Difference in energy')

figure('color',[1 1 1])
title('before')
yyaxis left
plot(t,(m_bal - m)./m);
% ylim([P_lower_bound-0.1 P_upper_bound+0.1])
xlabel('time [s]')
ylabel('Mass residual')
yyaxis right
plot(t,(E_bal - E)./E);
ylabel('Energy residual')
% ylim([v_min v_max])
% ylabel('v [m/s]')

%%
if strcmp(simType,'CAESCav')
    name = strcat(simType,'_P',num2str(P_L/1e6),'MPa_L',num2str(L),'m_Dt',num2str(floor(Dt/3600)),'h');
elseif strcmp(simType,'CAESPipe')
    name = strcat(simType,'_P',num2str(P_L/1e6),'MPa_L',num2str(L/1000),'km_Dt',num2str(floor(Dt/3600)),'h');
else
    disp('Unidentified simulation')
end

save(name)

%%
profile viewer
%%
% Pressure profile
figure('color',[1 1 1])
plot(x_n, P(:,2)./1e6)
hold all
plot(x_n, P(:,floor(n_t/5))./1e6)
plot(x_n, P(:,floor(2*n_t/5))./1e6)
plot(x_n, P(:,floor(3*n_t/5))./1e6)
plot(x_n, P(:,floor(4*n_t/5))./1e6)
plot(x_n, P(:,floor(n_t))./1e6)
ts = [t(2),t(floor(n_t/5)),t(floor(2*n_t/5)),t(floor(3*n_t/5)),t(floor(4*n_t/5)),t(n_t)]';
ts_leg = [num2str(ts),['s','s','s','s','s','s']'];
legend(ts_leg)
title('Pressure profiles [MPa]')

% Velocity profile
figure('color',[1 1 1])
plot(x_f, v(:,2))
hold all
plot(x_f, v(:,floor(n_t/5)))
plot(x_f, v(:,floor(2*n_t/5)))
plot(x_f, v(:,floor(3*n_t/5)))
plot(x_f, v(:,floor(4*n_t/5)))
plot(x_f, v(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Velocity profiles')

% density profile
figure('color',[1 1 1])
plot(x_n, rho(:,2))
hold all
plot(x_n, rho(:,floor(n_t/5)))
plot(x_n, rho(:,floor(2*n_t/5)))
plot(x_n, rho(:,floor(3*n_t/5)))
plot(x_n, rho(:,floor(4*n_t/5)))
plot(x_n, rho(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Density profiles')


%%
% Pressure field
figure('color',[1 1 1])
plot(t, P(2,:)./1e6)
hold all
plot(t, P(2*floor(n_n/5)+1,:)./1e6)
plot(t, P(3*floor(n_n/5)+1,:)./1e6)
plot(t, P(4*floor(n_n/5)+1,:)./1e6)
plot(t, P(n_n,:)./1e6)
xs = [x_n(2),x_n(2*floor(n_n/5)+1),x_n(3*floor(n_n/5)+1),x_n(4*floor(n_n/5)+1),x_n(n_n)]'
xs_leg = [num2str(xs),['m','m','m','m','m']'];
legend(xs_leg)
% legend('2*n/5','3*n/5','4*n/5','n')
title('Pressure x t')

% Pressure (faces) field
figure('color',[1 1 1])
plot(t, P_f(2,:)./1e6)
hold all
plot(t, P_f(2*floor(n_n/5)+1,:)./1e6)
plot(t, P_f(3*floor(n_n/5)+1,:)./1e6)
plot(t, P_f(4*floor(n_n/5)+1,:)./1e6)
plot(t, P_f(n_n,:)./1e6)
xs = [x_n(2),x_n(2*floor(n_n/5)+1),x_n(3*floor(n_n/5)+1),x_n(4*floor(n_n/5)+1),x_n(n_n)]'
xs_leg = [num2str(xs),['m','m','m','m','m']'];
legend(xs_leg)
% legend('2*n/5','3*n/5','4*n/5','n')
title('Pressure (faces) x t')

% Velocity field
figure('color',[1 1 1])
plot(t, v(2,:))
hold all
plot(t, v(2*floor(n_n/5),:))
plot(t, v(3*floor(n_n/5),:))
plot(t, v(4*floor(n_n/5),:))
plot(t, v(5*floor(n_n/5),:))
legend(xs_leg)
% legend('2*n/5','3*n/5','4*n/5','n')
title('Velocity x t')

% density (nodes) field
figure('color',[1 1 1])
plot(t, rho(2,:))
hold all
plot(t, rho(2*floor(n_n/5),:))
plot(t, rho(3*floor(n_n/5),:))
plot(t, rho(4*floor(n_n/5),:))
plot(t, rho(5*floor(n_n/5),:))
legend(xs_leg)
% legend('2*n/5','3*n/5','4*n/5','n')
title('Density x t')  

% density (faces) field
figure('color',[1 1 1])
plot(t, rho_f(2,:))
hold all
plot(t, rho_f(2*floor(n_n/5),:))
plot(t, rho_f(3*floor(n_n/5),:))
plot(t, rho_f(4*floor(n_n/5),:))
plot(t, rho_f(5*floor(n_n/5),:))
legend(xs_leg)
% legend('2*n/5','3*n/5','4*n/5','n')
title('Density (faces) x t')  

% Pressure profile
figure('color',[1 1 1])
plot(x_n, P(:,2)./1e6)
hold all
plot(x_n, P(:,floor(n_t/5))./1e6)
plot(x_n, P(:,floor(2*n_t/5))./1e6)
plot(x_n, P(:,floor(3*n_t/5))./1e6)
plot(x_n, P(:,floor(4*n_t/5))./1e6)
plot(x_n, P(:,floor(n_t))./1e6)
ts = [t(2),t(floor(n_t/5)),t(floor(2*n_t/5)),t(floor(3*n_t/5)),t(floor(4*n_t/5)),t(n_t)]';
ts_leg = [num2str(ts),['s','s','s','s','s','s']'];
legend(ts_leg)
title('Pressure profiles')

% Temperature profile
figure('color',[1 1 1])
plot(x_n, T(:,2))
hold all
plot(x_n, T(:,floor(n_t/5)))
plot(x_n, T(:,floor(2*n_t/5)))
plot(x_n, T(:,floor(3*n_t/5)))
plot(x_n, T(:,floor(4*n_t/5)))
plot(x_n, T(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Temperature profiles')

% Velocity profile
figure('color',[1 1 1])
plot(x_f, v(:,2))
hold all
plot(x_f, v(:,floor(n_t/5)))
plot(x_f, v(:,floor(2*n_t/5)))
plot(x_f, v(:,floor(3*n_t/5)))
plot(x_f, v(:,floor(4*n_t/5)))
plot(x_f, v(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Velocity profiles')

% density profile
figure('color',[1 1 1])
plot(x_n, rho(:,2))
hold all
plot(x_n, rho(:,floor(n_t/5)))
plot(x_n, rho(:,floor(2*n_t/5)))
plot(x_n, rho(:,floor(3*n_t/5)))
plot(x_n, rho(:,floor(4*n_t/5)))
plot(x_n, rho(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Density profiles')

% Pressure profile (faces)
figure('color',[1 1 1])
plot(x_f, P_f(:,2)./1e6)
hold all
plot(x_f, P_f(:,floor(n_t/5))./1e6)
plot(x_f, P_f(:,floor(2*n_t/5))./1e6)
plot(x_f, P_f(:,floor(3*n_t/5))./1e6)
plot(x_f, P_f(:,floor(4*n_t/5))./1e6)
plot(x_f, P_f(:,floor(n_t))./1e6)
ts = [t(2),t(floor(n_t/5)),t(floor(2*n_t/5)),t(floor(3*n_t/5)),t(floor(4*n_t/5)),t(n_t)]';
ts_leg = [num2str(ts),['s','s','s','s','s','s']'];
legend(ts_leg)
title('Pressure profiles (faces)')

% Temperature profile
figure('color',[1 1 1])
plot(x_f, T_f(:,2))
hold all
plot(x_f, T_f(:,floor(n_t/5)))
plot(x_f, T_f(:,floor(2*n_t/5)))
plot(x_f, T_f(:,floor(3*n_t/5)))
plot(x_f, T_f(:,floor(4*n_t/5)))
plot(x_f, T_f(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Temperature profiles (faces)')

% density profile
figure('color',[1 1 1])
plot(x_f, rho_f(:,2))
hold all
plot(x_f, rho_f(:,floor(n_t/5)))
plot(x_f, rho_f(:,floor(2*n_t/5)))
plot(x_f, rho_f(:,floor(3*n_t/5)))
plot(x_f, rho_f(:,floor(4*n_t/5)))
plot(x_f, rho_f(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Density profiles (faces)')

% Velocity profile
figure('color',[1 1 1])
plot(x_f, v(:,2))
hold all
plot(x_f, v(:,floor(n_t/5)))
plot(x_f, v(:,floor(2*n_t/5)))
plot(x_f, v(:,floor(3*n_t/5)))
plot(x_f, v(:,floor(4*n_t/5)))
plot(x_f, v(:,floor(n_t)))
legend(ts_leg)
% legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Velocity profiles (faces)')

if strcmp(simType,'CAESCav')
    name = strcat(simType,'_P',num2str(P_L/1e6),'MPa_L',num2str(L),'m_Dt',num2str(floor(Dt/3600)),'h');
elseif strcmp(simType,'CAESPipe')
    name = strcat(simType,'_P',num2str(P_L/1e6),'MPa_L',num2str(L/1000),'km_Dt',num2str(floor(Dt/3600)),'h');
else
    disp('Unidentified simulation')
end

save(name)
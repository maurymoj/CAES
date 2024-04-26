clear 

CP = py.importlib.import_module('CoolProp.CoolProp');
CP.PropsSI('D','P',101325,'T',298,'Air');

%--------------------- SIMULATION PARAMETERS ------------------------%
dx = 10;
dt = 1;

Dt = 10;

%----------------------- PROBLEM PARAMETERS -------------------------%
% Pipeline parameters
L = 50;
D = 0.5;
eps = 0.04e-3; % Absolute roughness 0.04 mm
epsD = eps/D;
A_h = pi*D^2/4;

% Inlet 
P_in = 101325;
T_in = 273.15 + 15;
m_dot = 1;

rho_in = CP.PropsSI('D','P',P_in,'T',T_in,'Air');
v_in = m_dot/(rho_in*A_h);

% Outlet
v_out = 0;

g = 9.81;
theta = 0;

%---------------------- ARRAYS INITIALIZATION ----------------------%
n = L/dx;       % number of nodes
N = n + 1;      % number of faces
n_t = Dt/dt+1;  % n of time steps

% Face initialization
v = zeros(N,n_t);
P_f = zeros(N,n_t);
T_f = zeros(N,n_t);
rho_f = zeros(N,n_t);

% Node initialization
v_n = zeros(n,n_t);
P = zeros(n,n_t);
T = zeros(n,n_t);
rho = zeros(n,n_t);

%------------------ SETTING INITIAL CONDITIONS ---------------------%
% Initial conditions at nodes (i -> x, j -> t)
% Thermodynamic properties over all solution space
P(:,1) = 101325*ones(n,1);
T(:,1) = 298*ones(n,1);
rho(:,1) = CP.PropsSI('D','P',P(1,1),'T',T(1,1),'Air');
% v_n = 0;

% Inlet Boundary condition
P_f(1,:) = P_in;
T_f(1,:) = T_in;
rho_f(1,:) = rho_in;

% Initial conditions at faces (i -> x, j -> t)
v(:,1) = 0; % V profile at t=0
v(:,1) = v_in; % V profile at t=0
v(1,:) = v_in; % Inlet flow rate (and velocity) is constant for all t>=0
v(end,:) = 0; % V is 0 at the pipe end for all t>=0

% Upwind scheme
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

% P_f(2:end,1) = P(:,1); % ASSUMING v >= 0 for t=0!!!!
% T_f(2:end,1) = T(:,1); % ASSUMING v >= 0 for t=0!!!!
% rho_f(2:end,1) = rho(:,1); % ASSUMING v >= 0 for t=0!!!!

T(:,:) = 298; % Isothermal pipeline assumption

P_corr = zeros(n,1);
rho_corr = zeros(n,1);
v_corr = zeros(N,1);






% ITERATION
j=2;

count = 0;
error_P = 10
while count < 100 && max(abs(P_corr./P(:,j))) > 0.001
    % Initial guess (properties at t+dt = properties at t)
    P(:,j) = P(:,j-1) + P_corr;
    rho(:,j) = rho(:,j-1) + rho_corr;
    T(:,j) = T(:,j-1);
    v(:,j) = v(:,j-1) + v_corr;
    
    



    % Upwind scheme
    v_n(:,j) = v_n(:,j-1);
    P_f(:,j) = P_f(:,j-1);
    T_f(:,j) = T_f(:,j-1);
    rho_f(:,j) = rho_f(:,j-1);
    




    %---------- v* calculation - Momentum control volume ---------------------%
    
    f = (2*log10(1/epsD)+1.14)^(-2); % Friction factor based on Nikuradse - IMPLEMENT COLEBROOK EQUATION
    
    % Matrix/vector initialization (only need to solve for the inner faces -
    % boundary faces are dealt with the boundary conditions)
    a = zeros(N-2);
    b = zeros(N-2,1);
    d = zeros(N-2,1);
    B = zeros(N-2,1);
    
    % For momentum equations the first line refers to face B and the last line to
    % face E, as we know velocities at faces A and F from the boundary
    % conditions
    a(1,2) = -max( 0, -rho(2,j)*v_n(2,j)/dx);    % a_C
    a(1,1) = rho_f(2,j)/dt + f*rho_f(2,j)*abs(v(2,j))/(2*D) ...
        + (rho(2,j)*v_n(2,j) - rho(1,j)*v_n(1,j))/dx ...
        - (a(1,2) - max(rho(1,j)*v_n(1,j)/dx,  0) );% a_B
    
    d(1) = -1/dx;
    b(1) = (rho_f(2,j-1)*v(2,j-1)/dt)-rho_f(2,j)*g*sind(theta)+max(rho(1,j)*v_n(1,j)/dx,  0)*v_n(1,j);
    
    B(1) = d(1)*(P(2,j)-P(1,j)) + b(1);
    
    
    a(end,end-1) = -max(rho(n-1,j)*v_n(n-1,j)/dx, 0); % a_N-2 (index N-2, N-3)
    % WHAT TO DO REGARDING TO THE a INDICES ??
    a(end,end) = ...                                  % a_N-1 (index N-2, N-2)
        rho_f(N-1,2)/dt + f*rho_f(N-1,j)*abs(v(N-1,j))/(2*D)...
        +(rho(n,j)*v_n(n,j) - rho(n-1,j)*v_n(n-1,j))/dx...
        - (a(end,end-1) + max(-rho(n,j)*v_n(n,j)/dx,  0) );
    
    d(end) = -1/dx;
    b(end) = (rho_f(N-1,j-1)*v(N-1,j-1)/dt)-rho_f(N-1,j)*g*sind(theta);
    
    B(end) = d(end)*(P(n,j)-P(n-1,j)) + b(end);
    
    for i=2:N-3
    % The indices of the node and face properties have to be added by 1, since
    % the 1st line was ommited (as we have it solved from the boundary
    % condition 
        a(i,i-1) = -max(rho(i,j)*v_n(i,j)/dx,                     0);
        a(i,i+1) = -max(                       0, -rho(i+1,j)*v_n(i+1,j)/dx);
        a(i,i)   = rho_f(i+1,j)/dt + f*rho_f(i+1,j)*abs(v(i+1,j))/(2*D) ...
            - (a(i,i-1) + a(i,i+1)) ...
            + (rho(i+1,j)*v_n(i+1,j)/dx - rho(i,j)*v_n(i,j)/dx);
        d(i) = -1/dx;
        b(i) = (rho_f(i+1,j-1)*v(i+1,j-1)/dt)-rho_f(i+1,j)*g*sind(theta);
        
        B(i) = d(i)*(P(i+1,j)-P(i,j)) + b(i);
    end
    
    % v_star = B\a
    v_star = linsolve(a,B);
    v_star = [v(1,j);v_star;v(end,j)]
    % v(:,j) = [v(1,j);v_star;v(end,j)];
    
    %--------------- Pressure correction P' calculations ---------------------%
    % The matrices/vectors in P' calculations have no shift, so index 1 means
    % control volume 1 and so on
    A = zeros(n);
    BB = zeros(n,1);
    u_sonic = zeros(N,1);
    u_sonic_n = zeros(n,1);
    drho_dP = zeros(N,1);
    drho_dP_n = zeros(n,1);
    
    for i = 1:n
        u_sonic(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'T',T_f(i,j),'Air');
        drho_dP(i) = 1/(u_sonic(i)^2);
    
        u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'T',T(i,j),'Air');
        drho_dP_n(i) = 1/(u_sonic_n(i)^2);
    end
    
    u_sonic(end) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'T',T_f(end,j),'Air');
    drho_dP(end) = 1/(u_sonic(end)^2);
    
    % left node
    % Coefficients of a and d are displaced by 1 (1 = B, 2 = C and so on) 
    A(1,2) = - max(0, -drho_dP(2)*v_star(2)/dx); % A_2
    A(1,1) = drho_dP_n(1)/dt - (rho_f(2,j)*d(1)/a(1,1))/dx + drho_dP(2)*v_star(2)/dx - A(1,2); % A_1
    
    BB(1) = (rho(1,j-1)-rho(1,j))/dt + (rho_f(1,j)*v_star(1) - rho_f(2,j)*v_star(2))/dx;
    
    
    % right node
    A(n,n-1) = (rho_f(N-1,j)*d(end)/a(end,end))/dx + max(drho_dP(N-1)*v_star(N-1)/dx, 0);
    A(n,n) = drho_dP_n(n)/dt - drho_dP(N-1)*v_star(N-1)/dx - A(n,n-1); % v_star(end) = 0
    
    BB(n) = (rho(n,j-1)-rho(n,j))/dt + rho_f(N-1,j)*v_star(N-1)/dx;
    
    for i=2:n-1
    % Subtracted 1 from a and d as there is no row for the first face -
    % boundary conditions
        A(i,i-1) = (rho_f(i,2)*d(i-1)/a(i-1,i-1))/dx - max(drho_dP(i)*v_star(i)/dx,                     0);
        A(i,i+1) = (rho_f(i+1,2)*d(i)/a(i,i))/dx - max(              0, -drho_dP(i+1)*v_star(i+1)/dx);
        A(i,i)   = drho_dP_n(i)/dt + (drho_dP(i+1)*v_star(i+1)/dx ...
            - drho_dP(i)*v_star(i)/dx) - ( A(i,i-1) + A(i,i+1) );
        
        BB(i) = (rho(i,j-1)-rho(i,j))/dt + (rho_f(i,j)*v_star(i) - rho_f(i+1,j)*v_star(i+1))/dx;
    end
    A
    BB
    P_corr = linsolve(A,BB)
    error = P_corr./P(:,j)
    
    %%
    rho_corr = drho_dP_n.*P_corr
    a_i = diag(a);
    v_corr = zeros(N,1);
    
    v_corr(1)         = d(1)./a_i(1).*P_corr(1);
    v_corr(2:end-1)   = d(2:end-1)./a_i(2:end-1).*(P_corr(2:end)-P_corr(1:end-1));
    v_corr(end)       = 0 % boundary condition
    % v_corr(end)       =  % Velocity corrected based on mass balance (m_in =
    % m_out)
    
    error_P = (P_corr./P(:,2))
    error_rho = (rho_corr./rho(:,2))
    error_v = (v_corr./v(:,2))

end


















































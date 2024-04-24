clear 

CP = py.importlib.import_module('CoolProp.CoolProp');
CP.PropsSI('D','P',101325,'T',298,'Air');

%--------------------- SIMULATION PARAMETERS ------------------------%
dx = 10;
dt = 1;

Dt = 10;

%----------------------- PROBLEM PARAMETERS -------------------------%
% Inlet 
P_in = 101325;
T_in = 273.15 + 15;
rho_in = CP.PropsSI('D','P',P_in,'T',T_in,'Air');
m_dot = 1;

% Outlet
v_out = 0;

% Pipeline properties
L = 50;
D = 0.5;
eps = 0.04e-3; % Absolute roughness 0.04 mm
epsD = eps/D;
A_h = pi*D^2/4;

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
P(:,1) = 101325*ones(n,1);
T(:,1) = 298*ones(n,1);
rho(:,1) = CP.PropsSI('D','P',P(1,1),'T',T(1,1),'Air');
% v_n = 0;

% Inlet Boundary condition
P_f(1,:) = P_in;
T_f(1,:) = T_in;
rho_f(1,:) = CP.PropsSI('D','P',P_f(1,1),'T',T_f(1,1),'Air');


% Upwind scheme
P_f(2:end,1) = P(:,1); % ASSUMING v >= 0 for t=0!!!!
T_f(2:end,1) = T(:,1); % ASSUMING v >= 0 for t=0!!!!
rho_f(2:end,1) = rho(:,1); % ASSUMING v >= 0 for t=0!!!!


% Initial conditions at faces (i -> x, j -> t)
v_in = m_dot/(rho_f(1,1)*A_h);
v(:,1) = 0;
% v(:,1) = v_in;
% v_n(:,1) = v_in;
% v(:,1) = v_in;
% v_n(:,1) = v_in;
v_n(:,1) = (v(1:end-1,1)>=0).*v((1:end-1),1) ...
    + (v(1:end-1,1)<0).*v((2:end),1);

% ITERATION

% Initial guess (properties at t+dt = properties at t
P(:,2) = P(:,1);
T(:,2) = T(:,1);
rho(:,2) = rho(:,1);
v_n(:,2) = v_n(:,1);

P_f(:,2) = P_f(:,1);
T_f(:,2) = T_f(:,1);
rho_f(:,2) = rho_f(:,1);
v(:,2) = v(:,1);

% v* calculation - Momentum control volume

f = (2*log10(1/epsD)+1.14)^(-2); % Friction factor based on Nikuradse - IMPLEMENT COLEBROOK EQUATION
a = zeros(N);
b = zeros(N,1);
d = zeros(N,1);
B = zeros(N,1);

a(1,2) = -max(0,-rho(1,2)*v_n(1,2)/dx);                                              % a_B
a(1,1) = rho_f(1,2)/dt + f*rho_f(1,2)*abs(v(1,2))/(2*D)+rho(1,2)*v_n(1,2)/dx-a(1,2); % a_A

d(1) = -1/dx;

b(1) = (rho_f(1,1)*v(1,1)/dt)+rho_in*v_in^2/dx-rho_f(1,2)*g*sind(theta);

B(1) = d(1)*(P(1,2)-P_in) + b(1);

a(end,end-1) = -max(rho(n,2)*v_n(n,2)/dx, 0); % a_B
a(end,end) = rho_f(end,2)/dt + f*rho_f(end,2)*abs(v(end,2))/(2*D)-rho(end,2)*v_n(end,2)/dx-a(end,end-1); % a_A

d(end) = -1/dx;
b(end) = (rho_f(end,1)*v(end,1)/dt)+rho_f(end,1)*v_out^2/dx-rho_f(end,2)*g*sind(theta);

B(end) = d(end)*(P(end,2)-P(end-1,2)) + b(end) % Should this line exist since we know v_end = 0 ?

for i=2:N-1
    a(i,i-1) = -max(rho(i-1,2)*v_n(i-1,2)/dx,                     0);
    a(i,i+1) = -max(                       0, -rho(i,2)*v_n(i,2)/dx);
    a(i,i)   = rho_f(i,2)/dt + f*rho_f(i,2)*abs(v(i,2))/(2*D) ...
        - (a(i,i-1) + a(i,i+1)) ...
        + (rho(i,2)*v_n(i,2)/dx - rho(i-1,2)*v_n(i-1,2)/dx);
    d(i) = -1/dx;
    b(i) = (rho_f(i,1)*v(i,1)/dt)-rho_f(i,2)*g*sind(theta);
    
    B(i) = d(i)*(P(i,2)-P(i-1,2)) + b(i);
end

% v_star = B\a
v_star = linsolve(a,B)


















































clear 
% close all
clc

CP = py.importlib.import_module('CoolProp.CoolProp');
CP.PropsSI('D','P',101325,'T',298,'Air');

%--------------------- SIMULATION PARAMETERS ------------------------%
dx = 100;
dt = 0.5;

Dt = 50;

%----------------------- PROBLEM PARAMETERS -------------------------%
% Pipeline parameters
L = 5000;
D = 0.75;
eps = 0.04e-3; % Absolute roughness 0.04 mm
epsD = eps/D;
A_h = pi*D^2/4;

% Initial conditions
P_o = 101325;
T_o = 273.15 + 15;

% A - Left side - Inlet 
P_in = 7000000;
T_in = 273.15 + 15;
m_in = 150;

rho_in = CP.PropsSI('D','P',P_in,'T',T_in,'Air');
cp_in = CP.PropsSI('C','P',P_in,'T',T_in,'Air');
v_in = m_in/(rho_in*A_h);


% B - Right side
% R_bound = 'Wall';
R_bound = 'Outlet';
% Outlet
if strcmp(R_bound,'Outlet')
    m_out = m_in;
elseif(strcmp(R_bound,'Wall'))
    m_out = 0;
    v_out = 0;
elseif(strcmp(R_bound,'Inlet'))
    
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

dm = m_in*dt;
dE = m_in*cp_in*T_in*dt;

g = 9.81;
theta = 0;

%---------------------- ARRAYS INITIALIZATION ----------------------%
n = L/dx+1;       % number of nodes
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

% Sanity check
m = zeros(n_t,1); % Total mass in pipeline
E = zeros(n_t,1); % Total energy in pipeline

%------------------ SETTING INITIAL CONDITIONS ---------------------%
% Initial conditions at nodes (i -> x, j -> t)
% Thermodynamic properties over all solution space
P(:,1) = P_o*ones(n,1);
T(:,1) = T_o*ones(n,1);
rho(:,1) = CP.PropsSI('D','P',P(1,1),'T',T(1,1),'Air');
cp = CP.PropsSI('C','P',P(1,1),'T',T(1,1),'Air');

% Inlet Boundary condition
P_f(1,:) = P_in;
T_f(1,:) = T_in;
rho_f(1,:) = rho_in;

% Initial conditions at faces (i -> x, j -> t)
v(:,1) = v_in; % V profile at t=0

v(1,:) = v_in; % Inlet flow rate (and velocity) is constant for all t>=0
% v(end,:) = v_out; % V is 0 at the pipe end for all t>=0

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

if strcmp(R_bound,'Outlet')
    v(end,1) = rho_f(1,1)/rho_f(end,1)*v(1,1); % ENSURING CONSERVATION OF MASS
elseif(strcmp(R_bound,'Wall'))
    v(end,1) = 0;
elseif(strcmp(R_bound,'Inlet'))
    
end 



T(:,:) = T_o; % Isothermal pipeline assumption

m(1) = sum(rho(:,1)*A_h*dx);
E(1) = sum(rho(:,1)*A_h*dx*cp.*T(:,1));

count = zeros(n_t,1);
f_guess = (2*log10(1/epsD)+1.14)^(-2); % Friction factor based on Nikuradse

for j=2:n_t
    % Initial guess for next time step is the same props as the previous t step
    P(:,j) = P(:,j-1);
    rho(:,j) = rho(:,j-1);
    T(:,j) = T(:,j-1);
    v(:,j) = v(:,j-1);
    
    P_corr = zeros(n,1);
    rho_corr = zeros(n,1);
    v_corr = zeros(N,1);
    
    v_star = v(:,j);
    count(j) = 0;
    error_P = 10;

    while count(j) < 100 && max(abs(error_P)) > 0.001

        alpha_P = 0.5;  % Pressure under-relaxation factor
        alpha_v = 0.5;  % velocity under-relaxation factor
        alpha_rho = 0.5;  % Density under-relaxation factor

        P(:,j) = P(:,j) + alpha_P*P_corr; % Pressure under-relaxation correction
        rho(:,j) = alpha_rho*(rho(:,j) + rho_corr) + (1-alpha_rho)*rho(:,j);

        v(1:end-1,j) = alpha_v*(v_star(1:end-1) + v_corr(1:end-1)) + (1-alpha_v)*v_star(1:end-1);
        
        if strcmp(R_bound,'Outlet')
            v(end,j) = rho(1,j)/rho(end,j)*v(1,j); % ENSURING CONSERVATION OF MASS
        elseif(strcmp(R_bound,'Wall'))
            v(end,j) = 0;
        % elseif(strcmp(R_bound,'Inlet'))
        end 

        % Initial guess (properties at t+dt = properties at t) - without
        % under-relaxation
        % P(:,j) = P(:,j) + P_corr;
        % rho(:,j) = rho(:,j) + rho_corr;
        % T(:,j) = T(:,j);
        % v(:,j) = v_star(:) + v_corr;
        
        % Upwind scheme
        v_n(:,j) = (v(1:end-1,j) >= 0).*v((1:end-1),j) ...
            +      (v(1:end-1,j) <  0).*v((2:end),j);
        
        % Properties in nodes to faces
        P_f(2:end-1,j) = (v(1:end-2,j) >= 0).*P(1:end-1,j) ...
            +            (v(1:end-2,j) <  0).*P(2:end,j);
        rho_f(2:end-1,j) = (v(1:end-2,j) >= 0).*rho(1:end-1,j) ...
            +            (v(1:end-2,j) <  0).*rho(2:end,j);

        P_f(end,j) = P(end,j); % ASSUMING v >= 0 for t>0!!!!
        rho_f(end,j) = rho(end,j); % ASSUMING v >= 0 for t>0!!!!


        % Properties
        u_sonic = zeros(N,1);
        drho_dP = zeros(N,1);
        nu = zeros(N,1);

        u_sonic_n = zeros(n,1);
        drho_dP_n = zeros(n,1);
        cp = zeros(n,1);

        % for i = 1:n % PROPERTIES FROM P AND T
        % 
        %     u_sonic(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'T',T_f(i,j),'Air');
        %     nu(i) = CP.PropsSI('viscosity','P',P_f(i,j),'T',T_f(i,j),'Air');
        %     cp(i) = CP.PropsSI('C','P',P_f(i,j),'T',T_f(i,j),'Air');
        % 
        %     drho_dP(i) = 1/(u_sonic(i)^2);
        % 
        %     u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'T',T(i,j),'Air');
        % 
        %     drho_dP_n(i) = 1/(u_sonic_n(i)^2);
        % end
        % 
        % u_sonic(end) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'T',T_f(end,j),'Air');
        % drho_dP(end) = 1/(u_sonic(end)^2);


        for i = 1:n  % PROPERTIES FROM P AND RHO
            T(i,j) = CP.PropsSI('T','P',P(i,j),'D',rho(i,j),'Air');
            cp(i) = CP.PropsSI('C','P',P(i,j),'D',rho(i,j),'Air');
            
            u_sonic(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),'Air');
            nu = CP.PropsSI('viscosity','P',P_f(i,j),'D',rho_f(i,j),'Air');

            drho_dP(i) = 1/(u_sonic(i)^2);

            u_sonic_n(i) = CP.PropsSI('speed_of_sound','P',P(i,j),'D',rho(i,j),'Air');
            drho_dP_n(i) = 1/(u_sonic_n(i)^2);
        end

        u_sonic(end) = CP.PropsSI('speed_of_sound','P',P_f(end,j),'D',rho_f(end,j),'Air');
        drho_dP(end) = 1/(u_sonic(end)^2);

        T_f(2:end-1,j) = (v(1:end-2,j) >= 0).*T(1:end-1,j) ...
            +            (v(1:end-2,j) <  0).*T(2:end,j);
        T_f(end,j) = T(end,j); % ASSUMING v >= 0 for t>0!!!!
    
        %---------- v* calculation - Momentum control volume ---------------------%
        
        % Friction factor based on Nikuradse - IMPLEMENT COLEBROOK EQUATION
        f = (2*log10(1/epsD)+1.14)^(-2)*ones(n,1);
        
        % Re = rho_f(:,j).*v(:,j)*D./nu(:);
        % f = zeros(n,1);
        % 
        % for i=1:n
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
        a = zeros(N-2);
        b = zeros(N-2,1);
        d = zeros(N-2,1);
        B = zeros(N-2,1);
        
        % For momentum equations the first line refers to face B and the last line to
        % face E, as we know velocities at faces A and F from the boundary
        % conditions
        a(1,2) = -max( 0, -rho(2,j)*v_n(2,j)/dx);    % a_C
        a(1,1) = rho_f(2,j)/dt + f(1)*rho_f(2,j)*abs(v(2,j))/(2*D) ...
            + (rho(2,j)*v_n(2,j) - rho(1,j)*v_n(1,j))/dx ...
            - (a(1,2) - max(rho(1,j)*v_n(1,j)/dx,  0) );% a_B
        
        d(1) = -1/dx;
        b(1) = (rho_f(2,j-1)*v(2,j-1)/dt)...
            - rho_f(2,j)*g*sind(theta)...
            + max(rho(1,j)*v_n(1,j)/dx,  0)*v_n(1,j);
        
        B(1) = d(1)*(P(2,j)-P(1,j)) ...
            + b(1);
        
        
        a(end,end-1) = -max(rho(n-1,j)*v_n(n-1,j)/dx, 0); % a_N-2 (index N-2, N-3)
        % WHAT TO DO REGARDING TO THE a INDICES ??
        a(end,end) = ...                                  % a_N-1 (index N-2, N-2)
            rho_f(N-1,j)/dt + f(N-1)*rho_f(N-1,j)*abs(v(N-1,j))/(2*D)...
            +(rho(n,j)*v_n(n,j) - rho(n-1,j)*v_n(n-1,j))/dx...
            - (a(end,end-1) - max(-rho(n,j)*v_n(n,j)/dx,  0) );
        
        d(end) = -1/dx;
        b(end) = (rho_f(N-1,j-1)*v(N-1,j-1)/dt)-rho_f(N-1,j)*g*sind(theta);
        
        B(end) = d(end)*(P(n,j)-P(n-1,j)) + b(end);
        
        for i=2:N-3
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
        % v(:,j) = [v(1,j);v_star;v(end,j)];

        %--------------- Pressure correction P' calculations ---------------------%
        % The matrices/vectors in P' calculations have no shift, so index 1 means
        % control volume 1 and so on
        A = zeros(n);
        BB = zeros(n,1);
        
        % left node
        % Coefficients of a and d are displaced by 1 (1 = B, 2 = C and so on) 
        A(1,2) = - max(0, -drho_dP(2)*v_star(2)/dx); % A_2
        A(1,1) = drho_dP_n(1)/dt - (rho_f(2,j)*d(1)/a(1,1))/dx + drho_dP(2)*v_star(2)/dx - A(1,2); % A_1
        
        BB(1) = (rho(1,j-1)-rho(1,j))/dt + (rho_f(1,j)*v_star(1) - rho_f(2,j)*v_star(2))/dx;
        
        
        % right node
        A(n,n-1) = (rho_f(N-1,j)*d(end)/a(end,end))/dx - max(drho_dP(N-1)*v_star(N-1)/dx, 0);
        A(n,n) = drho_dP_n(n)/dt - drho_dP(N-1)*v_star(N-1)/dx - A(n,n-1); % v_star(end) = 0
        
        BB(n) = (rho(n,j-1)-rho(n,j))/dt + (rho_f(N-1,j)*v_star(N-1) - rho_f(N,j)*v_star(N))/dx;
        
        for i=2:n-1
        % Subtracted 1 from a and d as there is no row for the first face -
        % boundary conditions for matrix/vector in momentum eq.
            A(i,i-1) = (rho_f(i,j)*d(i-1)/a(i-1,i-1))/dx - max(drho_dP(i)*v_star(i)/dx,                     0);
            A(i,i+1) = (rho_f(i+1,j)*d(i)/a(i,i))/dx - max(              0, -drho_dP(i+1)*v_star(i+1)/dx);
            A(i,i)   = drho_dP_n(i)/dt + (drho_dP(i+1)*v_star(i+1)/dx ...
                - drho_dP(i)*v_star(i)/dx) - ( A(i,i-1) + A(i,i+1) );
            
            BB(i) = (rho(i,j-1)-rho(i,j))/dt + (rho_f(i,j)*v_star(i) - rho_f(i+1,j)*v_star(i+1))/dx;
        end
        A;
        BB;
        P_corr = linsolve(A,BB);

        rho_corr = drho_dP_n.*P_corr;

        a_i = diag(a);
        v_corr = zeros(N,1);
        % v_corr(1)         = d(1)./a_i(1).*P_corr(1);
        v_corr(1)         = 0;% Inlet boundary condition
        v_corr(2:end-1)   = d./a_i.*(P_corr(2:end)-P_corr(1:end-1));
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        v_corr(end)       = 0; % Wall boundary condition
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        % v_corr(end)       =  % Velocity corrected based on mass balance (m_in =
        % m_out)

        error_P = (P_corr./P(:,j));
        error_rho = (rho_corr./rho(:,j));
        error_v = (v_corr./v(:,j));
        count(j) = count(j)+1;  

    end

    m(j) = sum(rho(:,j)*A_h*dx);
    E(j) = sum(rho(:,j)*A_h*dx.*cp.*T(:,j));
    [j count(j)]
end





mm = m(1) + dm*(0:n_t-1);
EE = E(1) + dE*(0:n_t-1);

x = 0:dx:L;
t = 0:dt:Dt;

% Figures of mass and energy over time
figure('color',[1 1 1]);plot(t,m)
hold on; plot(t,mm(1:end))
legend('m','\Delta m')
figure('color',[1 1 1]);plot(t,E)
hold on; plot(t,EE(1:end))
legend('E','\Delta E')

% differences between total mass/energy in and change in C.V. mass/energy
figure('color',[1 1 1]);plot(t,mm' - m)
title('Difference in mass')
figure('color',[1 1 1]);plot(t,EE' - E)
title('Difference in energy')
%

% Pressure field
figure('color',[1 1 1])
plot(t, P(2*floor(n/5),:))
hold all
plot(t, P(3*floor(n/5),:))
plot(t, P(4*floor(n/5),:))
plot(t, P(5*floor(n/5),:))
legend('2*n/5','3*n/5','4*n/5','n')
title('Pressure x t')

% Velocity field
figure('color',[1 1 1])
plot(t, v(2*floor(n/5),:))
hold all
plot(t, v(3*floor(n/5),:))
plot(t, v(4*floor(n/5),:))
plot(t, v(5*floor(n/5),:))
legend('2*n/5','3*n/5','4*n/5','n')
title('Velocity x t')

% density field
figure('color',[1 1 1])
plot(t, rho(2*floor(n/5),:))
hold all
plot(t, rho(3*floor(n/5),:))
plot(t, rho(4*floor(n/5),:))
plot(t, rho(5*floor(n/5),:))
legend('2*n/5','3*n/5','4*n/5','n')
title('Density x t')  

% Pressure profile
figure('color',[1 1 1])
plot(x, P(:,floor(n_t/5)))
hold all
plot(x, P(:,floor(2*n_t/5)))
plot(x, P(:,floor(3*n_t/5)))
plot(x, P(:,floor(4*n_t/5)))
plot(x, P(:,floor(n_t)))
legend('n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Pressure profiles')

% Temperature profile
figure('color',[1 1 1])
plot(x, T(:,2))
hold all
plot(x, T(:,floor(n_t/5)))
plot(x, T(:,floor(2*n_t/5)))
plot(x, T(:,floor(3*n_t/5)))
plot(x, T(:,floor(4*n_t/5)))
plot(x, T(:,floor(n_t)))
legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Temperature profiles')

% Velocity profile
figure('color',[1 1 1])
plot(x, v(1:end-1,2))
hold all
plot(x, v(1:end-1,floor(n_t/5)))
plot(x, v(1:end-1,floor(2*n_t/5)))
plot(x, v(1:end-1,floor(3*n_t/5)))
plot(x, v(1:end-1,floor(4*n_t/5)))
plot(x, v(1:end-1,floor(n_t)))
legend('dt','n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Velocity profiles')

% density profile
figure('color',[1 1 1])
plot(x, rho(:,floor(n_t/5)))
hold all
plot(x, rho(:,floor(2*n_t/5)))
plot(x, rho(:,floor(3*n_t/5)))
plot(x, rho(:,floor(4*n_t/5)))
plot(x, rho(:,floor(n_t)))
legend('n_t/5','2*n_t/5','3*n_t/5','4*n_t/5','n_t')
title('Density profiles')
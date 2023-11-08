% pressure over 70 km pipeline
pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter

%% Steady-state analysis with constant T - X km pipe - Panhandle B flow equation
clear
clc
% close all
%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15 + 273.15; % T = 15 oC - ~288 K
P_a = 101325 ; % P = 101.325 kPa
G = 1;      % Specific gas gravity - for air G = 1

P_in = 7e6; % 7 MPa
T_1 = 5 + 273.15; % Common operational condition assumption Nasr and Connor "Natural Gas Engineering and Safety Challenges"
T_f = T_1; % Isothermal assumption
    % !!! High temperatures on the outlet of compressor stations,
    % which can persist for up to 50 km. !!!

% Flow rate
Q_a = 14*1000000/(24*3600); % Conversion from mcmd (millions of cubic meters
    % per day to cubic meters per second) - real demand estimation 70 mscm/day
    % St Fergus, 14 is the proportional relative to area of one of the 3, 900
    % mm diameter pipes

% Pipe
D = 0.9; % 900 mm
dL = 1000;    % Distance increment
L_m = 70000; % Total distance [m]

eps = 0.04e-3; % Absolute roughness 0.04 mm

E = 0.75;   % Pipeline efficiency
    % Usually varies from 0.6 to 0.92 as a function of liquid presence and age
    % "Handbook of natural gas transmission and processing"

H_1 = 0;
H_2 = 130;

T_sin = 30000; % Period of sinusoidal height increment
dh_sin = 50;    % Amplitude of sinusoidal height increment

% Generate multiple plots while varying one variable var (can be any
% parameter set in the initial statements)
% Var = {1 50 100 200 300 400 500}; % Values to be taken by H
% Var = {130};                        % Value for H
%
Var = {15000 30000 60000}; % Values for height increment period
% Var = {0};
% Var = {0 25 50 75 100}; % Values for height increment amplitude

LStyle = {'b','r','k','b--','r--','k--','b-.'};

% % Sanity check with book equation comparisons
% G = 0.6
% T_a = (60 - 32)*5/9 + 273.15 % T = 15 oC - ~288 K
% P_a = 14.73*6894.76 % P = 14.73 psia
% rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
% L_m = 161000                           % m
% D = (16 - 2*0.250)*0.0254              % m
% Q_a = (100*0.02825)*1000000/(24*3600)     % Sm3/s
% T_f = (80 - 32)*5/9 + 273.15
% P_in = (1400 + 14.73)*6894.76
% eps = 0.7e-3*25.4e-3
% E = 0.95


%----------------------------------------------------------------------%

% Individual Figures
% P_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% U_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% psi_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% h_fig = figure('Color',[1 1 1]);
% hold on
% grid on
Psi_fig = figure('Color',[1 1 1]);
hold on
grid on


% 1 Figure with subplots
sp = figure('Color',[1 1 1]);
subplot(2,2,1)
hold on
grid on
subplot(2,2,2)
hold on
grid on
subplot(2,2,3)
hold on
grid on
subplot(2,2,4)
hold on
grid on

rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
T0 = T_a;
h0 = py.CoolProp.CoolProp.PropsSI('H','P',P_a,'T',T_a,'Air');
s0 = py.CoolProp.CoolProp.PropsSI('S','P',P_a,'T',T_a,'Air');

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
% gam = cp/cv;

load('StF_Ab2.mat');

for j=1:length(Var)
    % H_2 = Var{j};
    T_sin = Var{j};
    % dh_sin = Var{j};

    A = pi*D^2/4; % mm2
    L = 0:dL:L_m; % 70 km pipe with increments of 1 km

    % Elevation profile from St. Fergus -> Aberdeen
    H = interp1(x_SfA*1000,H_SfA,L,"linear","extrap");
    dh_SfA = diff(H);
    
    dh = (H_2-H_1)/length(L);

    D_mm = 1000*D;% mm
    epsD = eps/D;
    f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.
    
    m_dot = rho_a*Q_a;

    P_a_kPa = P_a/1000;
    Q_a_day = Q_a*24*3600;    % Sm3/day
    L_km = L./1000;           % km
    
    rho = zeros(length(L),1);
    nu= zeros(length(L),1);
    nu_Po = zeros(length(L),1);
    Z = zeros(length(L),1);
    u = zeros(length(L),1);
    Re = zeros(length(L),1);
    f = zeros(length(L),1);
    s = zeros(length(L),1);
    h = zeros(length(L),1);
    psi = zeros(length(L),1);
    U_erosional = zeros(length(L),1);
    P = zeros(length(L),1);
    T = T_f*ones(length(L),1);      % Assuming isothermal flow
    
    P(1) = P_in;
    
    for i=1:length(L)-1
        
        rho(i) = py.CoolProp.CoolProp.PropsSI('D','P',P(i),'T',T_f,'Air');
        nu(i) = py.CoolProp.CoolProp.PropsSI('V','P',P(i),'T',T_f,'Air');
        h(i) = py.CoolProp.CoolProp.PropsSI('H','P',P(i),'T',T_f,'Air');
        s(i) = py.CoolProp.CoolProp.PropsSI('S','P',P(i),'T',T_f,'Air');
        Z(i) = py.CoolProp.CoolProp.PropsSI('Z','P',P(i),'T',T_f,'Air');
        U_erosional(i) = 1.22*100/sqrt(rho(i)); % The constant 100, can be any value from 100 to 250, 

        nu_Po(i) = 10*nu(i);
    
        % Velocity m/s
        u(i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(i)*T_f/P(i)); % Q_b has to be converted to m3/day and D to mm
        if u(i) > 0.5*U_erosional(i)
            warning('U greater than 50% of the erosional velocity')
        end

        % Exergy
        % phi = h-h0 + T0*(s-s0)+u^2/2+gz
        % Negligible height difference
        psi(i) = h(i)-h0 - T0*(s(i)-s0)+u(i)^2/2;
        % Considering height difference
        % psi(i) = h(i)-h0 + T0*(s(i)-s0)+u(i)^2/2+g*z(i)
    
        Re(i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm
    
        % Friction factor
        % Iterations for the estimation of friction factor using Colebrook equation
        f_old = f_guess;
        df = 10;
        count=0;
        while (df > 0.0001 & count < 10) 
            % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
            f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
            df = (f_new - f_old)/f_old;
            f_old = f_new;
            count = count + 1; 
        end
        f(i) = f_old;
    
        P_1_kPa = P(i)/1000;
        P_2_kPa = 0.98*P_1_kPa; % Initial guess
        dP = 10;
        count = 0;
        
        while (dP > 0.0001 & count < 10)
            P_f_kPa = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
            Z_f = py.CoolProp.CoolProp.PropsSI('Z','P',P_f_kPa*1000,'T',T_f,'Air');
            % - Compressibility formula for natural gas
            % Z_f = 1/(1+
            % ((P_f*0.000145038-14.73)*344400*10^(1.785*G)/(T_f*1.8)^3.825));
            
            % dh = dh_SfA(i); % Custom elevation profile

            % P from General Flow Equation
            % Negligible height difference
            % P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(dL/1000)*Z_f*f(i) );
            % % Considering height difference
            % s_H = 0.0684*G*(dh)/(T_f*Z_f);                                % Constant gradient / Custom gradient 
            % s_H = 0.0684*G*(dh+dh_sin*sin(2*pi/T_sin*(i*dL) ))/(T_f*Z_f); % Sinusoidal profile
            % s_H = 0.0684*G*(dh))/(T_f*Z_f); % Profile St. F. -> Ab 
            % L_e = dL*(exp(s_H)-1)/s_H;
            % P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(L_e/1000)*Z_f*f(i) )/exp(s_H) );

            % P from Panhandle B equation
            % Negligible height difference
            % P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(dL/1000)*Z_f );
            % considering height difference
            % s_H = 0.0684*G*(dh)/(T_f*Z_f);                                % Constant gradient / Custom gradient
            s_H = 0.0684*G*(dh+dh_sin*sin(2*pi/T_sin*(i*dL) ))/(T_f*Z_f);   % Sinusoidal profile
            % s_H = 0.0684*G*(dh)/(T_f*Z_f);                                % Constant gradient / Custom gradient
            L_e = dL*(exp(s_H)-1)/s_H;
            % L_e = (dL/cos(atan(dh/dL)))*(exp(s_H)-1)/s_H;
            % L_e = (dL/cos( atan( (dh+dh_sin*sin(2*pi/T_sin*(i*dL)))/dL ))) *(exp(s_H)-1)/s_H;                   % Sinusoidal profile
            P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );
            % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
            if (Re(i) < 4e6) || (Re(i) > 40e6)
                warning('Re outside the indicated region for the Panhandle B equation.')
            end
            
            % CHECK CONDITIONS OF EQUATIONS

            dP = (P(i+1) - P_2_kPa)/P_2_kPa;
            P_2_kPa = P(i+1);
            count = count + 1;
        end
        P(i+1) = 1000*P(i+1); % conversion back to Pa
    
    end
    
    rho(end) = py.CoolProp.CoolProp.PropsSI('D','P',P(end),'T',T_f,'Air');
    h(end) = py.CoolProp.CoolProp.PropsSI('H','P',P(end),'T',T_f,'Air');
    s(end) = py.CoolProp.CoolProp.PropsSI('S','P',P(end),'T',T_f,'Air');
    nu(end) = py.CoolProp.CoolProp.PropsSI('V','P',P(end),'T',T_f,'Air');
    Z(end) = py.CoolProp.CoolProp.PropsSI('Z','P',P(end),'T',T_f,'Air');
    u(end) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(end)*T_f/P(end)); % Q_b has to be converted to m3/day and D to mm
    psi(end) = h(end)-h0 - T0*(s(end)-s0)+u(end)^2/2;
    
    % Plots
    
    % Individual plots
    % plot(L,P)           % Pa x m 
    % xlabel('L [m]')
    % ylabel('P [Pa]')
    % figure(P_fig)   % Pressure drop profile
    % plot(L/1000,P/1000) % kPa x km
    % plot(L/1000,P/1000,LStyle{j}) % kPa x km
    % xlabel('L [km]')
    % ylabel('P [kPa]')   
    % plot(L/1610,P*0.000145038 - 14.73) % psig x mi
    % xlabel('L [mi]')
    % ylabel('P [psig]')
    % ylim([1150 1450])
    % figure('Color',[1 1 1])
    % figure(U_fig)   % Velocity profile
    % plot(L/1000,u)
    % plot(L/1000,u,LStyle{j})
    % xlabel('L [km]')
    % ylabel('u [m/s]')
    % figure(psi_fig)   % Entropy profile
    % plot(L/1000,psi/1000,LStyle{j})
    % xlabel('L [km]')
    % ylabel('\Psi [kJ/kg]')
    % figure(h_fig) % Height profile
    % plot(L/1000, (1:length(L))*dh + dh_sin*sin(2*pi/T_sin*dL*(1:length(L))))
    % xlabel('L [km]')
    % ylabel('H [m]')
    figure(Psi_fig)
    plot(L/1000,m_dot*psi/1000,LStyle{j})
    xlabel('L [km]')
    ylabel('\Psi [kW]')

    % Subplot structure
    figure(sp)
    subplot(2,2,1)
    plot(L/1000,P/1000,LStyle{j})
    xlabel('L [km]')
    ylabel('P [kPa]') 
    % Velocity profile
    subplot(2,2,2)
    plot(L/1000,u,LStyle{j})
    xlabel('L [km]')
    ylabel('u [m/s]')
    subplot(2,2,3)
    % Entropy profile
    plot(L/1000,psi/1000,LStyle{j})
    xlabel('L [km]')
    ylabel('\psi [kJ/kg]')
    subplot(2,2,4)
    % Height profile
    % plot(L/1000, (1:length(L))*dh,LStyle{j})
    % plot(L/1000, H,LStyle{j})
    plot(L/1000, (1:length(L))*dh + dh_sin*sin(2*pi/T_sin*dL*(1:length(L))),LStyle{j})
    % plot(L/1000, H,LStyle{j})
    xlabel('L [km]')
    ylabel('H [m]')

end

% figure(P_fig)
% legend(string(Var))
% 
% figure(U_fig)
% legend(string(Var))
% 
% figure(psi_fig)
% legend(string(Var))
% 
% figure(h_fig)
% legend(string(Var))

figure(sp)
legend(string(Var),'Location','northwest')
% figure('Color',[1 1 1]) % Height profile
% plot(L/1000, (1:length(L))*dh + dh_sin*sin(2*pi/30000*dL*(1:length(L))))
% xlabel('L [km]')
% ylabel('H [m]')



%% Steady-state analysis with variable T - 1 km pipe - Panhandle B flow equation
clear
clc
% close all
%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15+273.15; % T = 15 oC - ~288 K
P_a = 101325 ; % P = 101.325 kPa
G = 1;      % Specific gas gravity - for air G = 1

P_in = 7e6; % 7 MPa
T_in = 5 + 273.15; % Compressor output temperature

% Flow rate
Q_a = 14*1000000/(24*3600); % Conversion from mcmd (millions of cubic meters
    % per day to cubic meters per second) - real demand estimation 70 mscm/day
    % St Fergus, 14 is the proportional relative to area of one of the 3, 900
    % mm diameter pipes

% Pipe
D = 0.9; % 900 mm
dL = 1000;    % Distance increment [m]
L_m = 70000; % Total distance [m]
% L_m = 240000; % Total distance [m]

eps = 0.04e-3; % Absolute roughness 0.04 mm

E = 0.75;   % Pipeline efficiency
    % Usually varies from 0.6 to 0.92 as a function of liquid presence and age
    % "Handbook of natural gas transmission and processing"

% Elevation
H_1 = 0;
H_2 = 130;

% Soil temperature, assumed constant along the pipeline
T_s = 5 + 273.15; % 5 °C (Nasr and Connor "Natural Gas Engineering and Safety Challenges")


% Generate multiple plots while varying one variable var (can be any
% parameter set in the initial statements)
% Var = {10 20 50 100 200 500}; % Values to be taken by var - Height
% Var = {278.15 283 288 293 298 323 373}; % Values to be taken by var - Inlet Temperature
% Var = {278.15 298 318 338 358 378}; % Values to be taken by var - Inlet Temperature
% Var = {100 + 273.15};
Var = {50 + 273.15};
LStyle = {'b','r','k','b--','r--','k--','b-.'};

%----------------------------------------------------------------------%

% P_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% U_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% psi_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% psi_comp_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% h_fig = figure('Color',[1 1 1]);
% hold on
% grid on
% T_fig = figure('Color',[1 1 1]);
% hold on
% grid on

% 1 Figure with subplots
sp = figure('Color',[1 1 1]);
subplot(2,2,1)
hold on
grid on
subplot(2,2,2)
hold on
grid on
subplot(2,2,3)
hold on
grid on
subplot(2,2,4)
hold on
grid on

rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
T0 = T_a;
h0 = py.CoolProp.CoolProp.PropsSI('H','P',P_a,'T',T_a,'Air');
s0 = py.CoolProp.CoolProp.PropsSI('S','P',P_a,'T',T_a,'Air');

cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;

m_dot = rho_a*Q_a;
% j=1;
load('StF_Ab2.mat');

for j=1:length(Var)
    % H_2 = Var{j};
    T_in = Var{j};

    A = pi*D^2/4; % mm2
    L = 0:dL:L_m; % 70 km pipe with increments of 1 km
    H = interp1(x_SfA*1000,H_SfA,L,"linear","extrap");
    dh_SfA = diff(H);    
    
    dh = (H_2-H_1)/length(L);
    D_mm = 1000*D;% mm
    
    epsD = eps/D;
    f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation using Nikuradse eq.
    
    m_dot = rho_a*Q_a;

    P_a_kPa = P_a/1000;
    Q_a_day = Q_a*24*3600;    % Sm3/day
    L_km = L./1000;           % km
    
    rho = zeros(length(L),1);
    nu= zeros(length(L),1);
    nu_Po = zeros(length(L),1);
    Z = zeros(length(L),1);
    u = zeros(length(L),1);
    Re = zeros(length(L),1);
    f = zeros(length(L),1);
    s = zeros(length(L),1);
    h = zeros(length(L),1);
    psi = zeros(length(L),1);
    U_erosional = zeros(length(L),1);
    P = zeros(length(L),1);
    T = T_a*ones(length(L),1);
    
    
    P(1) = P_in;
    T(1) = T_in;
    
    for i=1:length(L)-1
        dT = 10;
        count_T = 0;
        T_old = max(0.98*T(i),T_s);
        while (dT > 0.0001 & count_T < 1000) 
            
            if T_old == T_s
                T_f = T_s;
            else
                T_f = T_s + (T(i)-T_old)/log((T(i)-T_s)/(T_old-T_s)); % assuming soil T below gas temp and joule-thompson effect negligible
            end

            rho(i) = py.CoolProp.CoolProp.PropsSI('D','P',P(i),'T',T_f,'Air');
            nu(i) = py.CoolProp.CoolProp.PropsSI('V','P',P(i),'T',T_f,'Air');
            h(i) = py.CoolProp.CoolProp.PropsSI('H','P',P(i),'T',T_f,'Air');
            s(i) = py.CoolProp.CoolProp.PropsSI('S','P',P(i),'T',T_f,'Air');
            Z(i) = py.CoolProp.CoolProp.PropsSI('Z','P',P(i),'T',T_f,'Air');
            U_erosional(i) = 1.22*100/sqrt(rho(i)); % The constant 100, can be any value from 100 to 250, 

            nu_Po(i) = 10*nu(i);
        
            % Velocity m/s
            u(i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(i)*T_f/P(i)); % Q_b has to be converted to m3/day and D to mm
            if u(i) > 0.5*U_erosional(i)
                warning('U greater than 50% of the erosional velocity')
            end

            % Exergy
            % phi = h-h0 + T0*(s-s0)+u^2/2+gz
            % Negligible height difference
            psi(i) = h(i)-h0 - T0*(s(i)-s0)+u(i)^2/2;
            % Considering height difference
            % phi(i) = h(i)-h0 + T0*(s(i)-s0)+u(i)^2/2+g*z(i)
        
            Re(i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P converted to kPa, Q to m3/day, nu to poise (1 Pa s = 10 poise), and D to mm
        
            % Friction factor
            % Iterations for the estimation of friction factor using Colebrook equation
            f_old = f_guess;
            df = 10;
            count=0;
            while (df > 0.0001 & count < 10) 
                % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
                f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
                df = (f_new - f_old)/f_old;
                f_old = f_new;
                count = count + 1; 
            end
            f(i) = f_old;
        
            P_1_kPa = P(i)/1000;
            P_2_kPa = 0.98*P_1_kPa; % Initial guess
            dP = 10;
            count = 0;
            
            while (dP > 0.0001 & count < 1000)
                P_f_kPa = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
                % P_f = (P_1_kPa + P_2_kPa)/2;
                Z_f = py.CoolProp.CoolProp.PropsSI('Z','P',P_f_kPa*1000,'T',T_f,'Air');
                % Z_f = 1/(1+
                % ((P_f*0.000145038-14.73)*344400*10^(1.785*G)/(T_f*1.8)^3.825));
                % - Compressibility formula for natural gas

                dh = dh_SfA(i); % Custom elevation profile

                % P from General Flow Equation
                % Negligible height difference
                % P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(dL/1000)*Z_f*f(i) );
                % % Considering height difference
                % s_H = 0.0684*G*(dh)/(T_f*Z_f);                                % Constant gradient / Custom gradient 
                % s_H = 0.0684*G*(dh+dh_sin*sin(2*pi/T_sin*(i*dL) ))/(T_f*Z_f); % Sinusoidal profile
                % s_H = 0.0684*G*(dh))/(T_f*Z_f); % Profile St. F. -> Ab 
                % L_e = dL*(exp(s_H)-1)/s_H;
                % P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(L_e/1000)*Z_f*f(i) )/exp(s_H) );

                % P from Panhandle B equation
                % Negligible height difference
                % P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(dL/1000)*Z_f );
                % considering height difference
                s_H = 0.0684*G*(dh)/(T_f*Z_f);
                L_e = dL*(exp(s_H)-1)/s_H;
                P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );
                % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                if (Re(i) < 4e6) || (Re(i) > 40e6)
                    warning('Re outside the indicated region for the Panhandle B equation.')
                end

                dP = (P(i+1) - P_2_kPa)/P_2_kPa;
                P_2_kPa = P(i+1);
                count = count + 1;
            end
            P(i+1) = 1000*P(i+1); % conversion back to Pa
            
            %-----------------------------------------------------------------
            % Temperature profile
            % Convection heat transfer assuming constant pipe temperature (5 °C)
            Pr = py.CoolProp.CoolProp.PropsSI('Prandtl','P',P_f_kPa*1000,'T',T_f,'Air');
            k_fl = py.CoolProp.CoolProp.PropsSI('conductivity','P',P_f_kPa*1000,'T',T_f,'Air');

            % Bhatti and Shah
            % Appropriate for rough pipes and Re > 10000, 0.5 < Pr < 10,
            % 0.002 < e/D < 0.05
            Re_e=Re(i)*epsD*sqrt(f(i)/8);
            Nu = ( (f(i)/8)*Re(i)*Pr )/( 1+sqrt(f(i)/8)*(4.5*Re_e^0.2*sqrt(Pr)-8.48) );
            if Re(i) < 10000
                warning('Re < 10000, Nu correlation not appropriate')
            elseif (Pr < 0.5) || (Pr > 10)
                warning('Pr outside the Nu correlation appropriate region (0.5 < Pr < 10)')
            end

            % Gnielinski Nu equation, 
            % Appropriate for 0.5 < Pr < 2000, 3000 < Re < 5e6
            % Adequate for a first approx. for rough tubes, but not ideal
            % Nu = ( (f(i)/8)*(Re(i)-1000)*Pr )/( 1+12.7*(f(i)/8)^0.5*(Pr^(2/3)-1) );
            % if (Re(i) < 3000) || (Re(i) > 5e6)
            %     warning('Re outside the Nu correlation appropriate region (3000 < Re < 5e6)')
            % elseif (Pr < 0.5) || (Pr > 2000)
            %     warning('Pr outside the Nu correlation appropriate region (0.5 < Pr < 2000)')
            % end
            
            % Adjusted Colburn Nu equation for cooling - "Themal fluid sciences - Cengel"
            % Appropriate for 0.6 < Pr < 160, Re > 10000
            % Nu = 0.023*Re(i)^0.8*Pr^0.4; % Adjusted Colburn equation for cooling - "Thermal fluid sciences - Cengel" not ideal for non-smooth pipes!!
            % if (Re(i) < 10000)
            %     warning('Re outside the Nu correlation appropriate region (Re > 10000)')
            % elseif (Pr < 0.6) || (Pr > 160)
            %     warning('Pr outside the Nu correlation appropriate region (0.6 < Pr < 160)')
            % end

            % Nu_C(i) = Nu;
            % Nu = Nu_B(i);
            U = Nu*k_fl/D;
            
            %-----------------------------------------------------------------
            C_p = py.CoolProp.CoolProp.PropsSI('C','P',P_f_kPa*1000,'T',T_f,'Air');
            a = pi*D*U/(m_dot*C_p);
            T(i+1) = T_s + (T(i) - T_s)*exp(-a*dL);
            dT = (T(i+1) - T_old)/T_old;
            T_old = max(T(i+1),T_s);
            count_T = count + 1; 
        end
    end
    
    if T(end) == T_s
       T_f = T_s;
    else
       T_f = T_s + ( T(end-1)-T(end) )/log( (T(end-1)-T_s)/(T(end)-T_s) ); % assuming soil T below gas temp and joule-thompson effect negligible
    end

    rho(end) = py.CoolProp.CoolProp.PropsSI('D','P',P(end),'T',T_f,'Air');
    h(end) = py.CoolProp.CoolProp.PropsSI('H','P',P(end),'T',T_f,'Air');
    s(end) = py.CoolProp.CoolProp.PropsSI('S','P',P(end),'T',T_f,'Air');
    nu(end) = py.CoolProp.CoolProp.PropsSI('V','P',P(end),'T',T_f,'Air');
    Z(end) = py.CoolProp.CoolProp.PropsSI('Z','P',P(end),'T',T_f,'Air');
    u(end) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(end)*T_f/P(end)); % Q_b has to be converted to m3/day and D to mm
    psi(end) = h(end)-h0 - T0*(s(end)-s0)+u(end)^2/2;
    
    % Plots
    % figure('Color',[1 1 1])
    % plot(L,P)           % Pa x m 
    % xlabel('L [m]')
    % ylabel('P [Pa]')
    % figure(P_fig)   % Pressure drop profile
    % plot(L/1000,P/1000,LStyle{j}) % kPa x km
    % xlabel('L [km]')
    % ylabel('P [kPa]')
    % plot(L/1610,P*0.000145038 - 14.73) % psig x mi
    % xlabel('L [mi]')
    % ylabel('P [psig]')
    % ylim([1150 1450])
    
    % figure('Color',[1 1 1])
    % figure(U_fig)   % Velocity profile
    % plot(L/1000,u,LStyle{j})
    % xlabel('L [km]')
    % ylabel('u [m/s]')
    
    % figure(psi_fig) % Exergy profile
    % plot(L/1000,psi/1000,LStyle{j})
    % % plot(L/1000,phi,LStyle{j})
    % xlabel('L [km]')
    % ylabel('\Psi [kJ/kg]')

    % figure(T_fig)   % Temperature profile
    % plot(L/1000,T,LStyle{j})
    % xlabel('L [km]')
    % ylabel('T [K]')
    % 
    % figure(psi_comp_fig) % Exergy profile
    % plot(L/1000,(h-h0)./1000,LStyle{1})
    % plot(L/1000,(-T0*(s-s0))./1000,LStyle{2})
    % plot(L/1000,(u.^2/2)./1000,LStyle{3})
    % title('T_{in} = ',T_in)
    % legend('h-h0','-T0 (s-s0)','u^2/2')
    % 
    % % plot(L/1000,phi,LStyle{j})
    % xlabel('L [km]')
    % ylabel('\Psi components[kJ/kg]')

    % Subplot structure
    figure(sp)
    subplot(2,2,1)
    plot(L/1000,P/1000,LStyle{j}) % kPa x km
    xlabel('L [km]')
    ylabel('P [kPa]') 
    % Velocity profile
    subplot(2,2,2)
    plot(L/1000,u,LStyle{j})
    xlabel('L [km]')
    ylabel('u [m/s]')
    subplot(2,2,3)
    % Entropy profile
    plot(L/1000,psi/1000,LStyle{j})
    xlabel('L [km]')
    ylabel('\Psi [kJ/kg]')
    subplot(2,2,4)
    % Height profile
    % plot(L/1000, (1:length(L))*dh,LStyle{j})
    % plot(L/1000, (1:length(L))*dh + dh_sin*sin(2*pi/T_sin*dL*(1:length(L))),LStyle{j})
    plot(L/1000, H,LStyle{j})
    xlabel('L [km]')
    ylabel('H [m]')
end

% figure(P_fig)
% legend(string(Var))
% 
% figure(U_fig)
% legend(string(Var))
% 
% figure(psi_fig)
% legend(string(Var))
% 
% figure(T_fig)
% legend(string(Var))
% 
% figure(psi_comp_fig)
% legend('h-h0','-T0 (s-s0)','u^2/2')

figure(sp)
legend(string(Var),'Location','northwest')
%% Erosional velocity Ve
% usually velocity should be limited to 20 m/s. Operational velocity is
% usually limited to 50% of Ve
% Ambient conditions
T_a = 15+273.15; % T = 15 oC - ~288 K
P_a = 101325 ; % P = 101.325 kPa
T = 323;
P = 7000000;

rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
rho = py.CoolProp.CoolProp.PropsSI('D','P',P,'T',T,'Air');

cte = 100;  % cte ranges from 100 to 250
            % cte = 100 for continuous service
            %       125 for noncontinuous service
            %       120-200 for continuous, noncorrosive or corrosive
            %       controlled, if no solid particles are present
            % "Handbook of Natural Gas Transmission and Processing"
Ve_a = 1.22*cte/sqrt(rho_a) 
Ve = 1.22*cte/sqrt(rho) 
% Equation in Nasr and Connor, 2014 - Natural Gas Engineering and Safety 
% Challenges

%% Gas temperature profile
LStyle = ( 2*pi*(D/2)*U )/( m_dot*c_p )
T_X = T_s + (T_1-T_s)*exp(-LStyle*X)
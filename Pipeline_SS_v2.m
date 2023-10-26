pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter

%% Steady-state analysis with constant T - 1 km pipe - Panhandle B flow equation
clear
clc
% close all
%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15 + 273.15;  % T = 15 oC - ~288 K
P_a = 101325 ;      % P = 101.325 kPa
G = 1;              % Specific gas gravity - for air G = 1

P_in = 7e6;         % 7 MPa
T_1 = 5 + 273.15;   % Common operational condition assumption Nasr and Connor "Natural Gas Engineering and Safety Challenges"
T_f = T_1;          % Isothermal assumption
              % !!! High temperatures on the outlet of compressor stations,
              % which can persist for up to 50 km. !!!

% Flow rate
Q_a = 14*1000000/(24*3600); % Conversion from mcmd (millions of cubic meters
    % per day to cubic meters per second) - real demand estimation 70 mscm/day
    % St Fergus, 14 is the proportional relative to area of one of the 3, 900
    % mm diameter pipes

% Pipe
D = 0.9;      % 900 mm
dL = 1000;    % Distance increment [m]
L_m = 120000; % Total distance [m]

eps = 0.04e-3; % Absolute roughness 0.04 mm - excel

E = 0.75;   % Pipeline efficiency
    % Usually varies from 0.6 to 0.92 as a function of liquid presence and age
    % "Handbook of natural gas transmission and processing"

H_1 = 0; % elevation at pipe inlet
H_2 = 1; % elevation at pipe outlet

% Generate multiple plots while varying one variable var (can be any
% parameter set in the initial statements)
% Var = {1 50 100 200 300 400 500}; % Values to be taken by var H2
Var = {1 50 100 150 200 250 300 500}; % Values to be taken by var H2
% Var = {1};

% Plot styles
LStyle = {'b','r','k','b--','r--','k--','b-.','r-.','k-.'};

%----------------------------------------------------------------------%

% Creation of figures for the Pressure, Velocity (U) and Exergy plots
P_fig = figure('Color',[1 1 1]);
hold on
grid on
U_fig = figure('Color',[1 1 1]);
hold on
grid on
psi_fig = figure('Color',[1 1 1]);
hold on
grid on

% Referece state parameters
T0 = T_a;
h0 = py.CoolProp.CoolProp.PropsSI('H','P',P_a,'T',T_a,'Air');
s0 = py.CoolProp.CoolProp.PropsSI('S','P',P_a,'T',T_a,'Air');
rho_a = py.CoolProp.CoolProp.PropsSI('D','P',P_a,'T',T_a,'Air');
cp = py.CoolProp.CoolProp.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = py.CoolProp.CoolProp.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
gam = cp/cv;


for j=1:length(Var)
    H_2 = Var{j};

    A = pi*D^2/4;               % Pipe diameter [m2]
    L = 0:dL:L_m;               % Discretization of pipe in increments of dL km
    dh = (H_2-H_1)/length(L);   % Constant Height increment
    D_mm = 1000*D;              % Pipe diameter [mm]
    
    epsD = eps/D;               % Relative roughness
    f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial estimation of f using Nikuradse eq.
    
    P_a_kPa = P_a/1000;         % Ambient pressure [kPa]
    Q_a_day = Q_a*24*3600;      % Volumetric flow rate [Sm3/day]
    L_km = L./1000;             % Length [km]
    
    % Vector initialization
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
    P = zeros(length(L),1);

    T = T_f*ones(length(L),1);  % Isothermal flow assumption
    
    P(1) = P_in;
    
    for i=1:length(L)-1
        
        rho(i) = py.CoolProp.CoolProp.PropsSI('D','P',P(i),'T',T_f,'Air');
        nu(i) = py.CoolProp.CoolProp.PropsSI('V','P',P(i),'T',T_f,'Air');
        h(i) = py.CoolProp.CoolProp.PropsSI('H','P',P(i),'T',T_f,'Air');
        s(i) = py.CoolProp.CoolProp.PropsSI('S','P',P(i),'T',T_f,'Air');
        Z(i) = py.CoolProp.CoolProp.PropsSI('Z','P',P(i),'T',T_f,'Air');
        
        nu_Po(i) = 10*nu(i);
    
        % Velocity m/s
        u(i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(i)*T_f/P(i)); % Q_b has to be converted to m3/day and D to mm
        
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
        
        while (dP > 0.0001 & count < 10)
            P_f_kPa = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
            % P_f = (P_1_kPa + P_2_kPa)/2;
            Z_f = py.CoolProp.CoolProp.PropsSI('Z','P',P_f_kPa*1000,'T',T_f,'Air');
            % Z_f = 1/(1+
            % ((P_f*0.000145038-14.73)*344400*10^(1.785*G)/(T_f*1.8)^3.825));
            % - Compressibility formula for natural gas
        
            % P from Panhandle B equation
            % Negligible height difference
            % P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(dL/1000)*Z_f );
            % considering height difference
            s_H = 0.0684*G*(dh)/(T_f*Z_f);
            L_e = dL*(exp(s_H)-1)/s_H;
            P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );
            
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
    psi(end) = h(end)-h0 - T0*(s(end)-s0)+u(end)^2/2; % +gz ?
    
    % Plots
    
    % plot(L,P)           % Pa x m 
    % xlabel('L [m]')
    % ylabel('P [Pa]')
    figure(P_fig)   % Pressure drop profile
    plot(L/1000,P/1000,LStyle{j}) % kPa x km
    xlabel('L [km]')
    ylabel('P [kPa]')   
    % plot(L/1610,P*0.000145038 - 14.73) % psig x mi
    % xlabel('L [mi]')
    % ylabel('P [psig]')
    % ylim([1150 1450])
    
    % figure('Color',[1 1 1])
    figure(U_fig)   % Velocity profile
    plot(L/1000,u,LStyle{j})
    xlabel('L [km]')
    ylabel('u [m/s]')
    
    figure(psi_fig)   % Entropy profile
    plot(L/1000,psi/1000,LStyle{j})
    xlabel('L [km]')
    ylabel('\Psi [kJ/kg]')
end

figure(P_fig)
legend(string(Var))

figure(U_fig)
legend(string(Var))

figure(psi_fig)
legend(string(Var))
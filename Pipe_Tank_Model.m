% pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter

%%
clear
clc
close all
CP = py.importlib.import_module('CoolProp.CoolProp');

%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15 + 273.15; % T = 15 oC - ~288 K
P_a = 101325 ; % P = 101.325 kPa
G = 1;      % Specific gas gravity - for air G = 1

P_in = 7e6; % 7 MPa
T_1 = 5 + 273.15; % Common operational condition assumption Nasr and Connor "Natural Gas Engineering and Safety Challenges"
T_f = T_1; % Isothermal assumption

% Pipe
D = 0.9; % 900 mm
% D = 1.2;
% D = 0.45;
dL = 1000;    % Distance increment
L_m = 400000; % Total distance [m]

eps = 0.04e-3; % Absolute roughness 0.04 mm

E = 0.75;   % Pipeline efficiency
    % Usually varies from 0.6 to 0.92 as a function of liquid presence and age
    % "Handbook of natural gas transmission and processing"

H_1 = 0;
H_2 = 1;

Vol_tank = pi*D^2/4 * L_m; % Pipeline volume

% Elevation profile
% H_prof = "Horizontal";
H_prof = "Fixed tilt";
% H_prof = "Sinusoidal";
% H_prof = "Custom_prof";

% Flow equation
% Flow_eq = "GFE";
Flow_eq = "PanB";

% Transient parameters
dt = 10; % s
T_total = 3600;
% t = 0:dt:T_total;
t = 0;

% Generate multiple plots while varying one variable var (can be any
% parameter set in the initial statements)
% Pipe diameter mm
% Var = {0.450,0.900,1.200};
% Var = {0.900,1.200};
% Flow rate
% Var2 = {3*1000000/(24*3600),14*1000000/(24*3600),25*1000000/(24*3600)};
% Elevation at H2
Var = {100};
% Var = {1 50 100 200 300 400 500}; % Values to be taken by H
% Var = {50 100 150 200 250 300};     % Value for H
% Sinusoidal period - height increment
% Var = {5000 10000 25000 50000};
% Sinusoidal amplitude - height increment
% Var = {0 25 50 75 100}; % Values for height increment amplitude
% Elevation profile
% Var = {"Horizontal","Fixed tilt","Sinusoidal","Custom prof"};
% Var = {"Horizontal","Fixed tilt","Custom prof"};
% Flow equation
% Var = {"GFE","PanB"};

LStyle = {'b','r','k','b--','r--','k--','b-.'};

%----------------------------------------------------------------------%

% Individual Figures
P_fig = figure('Color',[1 1 1]);
hold on
grid on
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
pct_Psi_fig = figure('Color',[1 1 1]);
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

rho_a = CP.PropsSI('D','P',P_a,'T',T_a,'Air');
T0 = T_a;
h0 = CP.PropsSI('H','P',P_a,'T',T_a,'Air');
s0 = CP.PropsSI('S','P',P_a,'T',T_a,'Air');

cp = CP.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = CP.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
% gam = cp/cv;

% load('StF_Ab2.mat');

for k = 1:length(Var)
    
    A = pi*D^2/4; % Area, mm2
    % L = 0:dL:L_m; % pipe length with increments of dL, m
    % L = [0;dL]; % pipe length with increments of dL, m
    L = 0;

    % Elevation profile from St. Fergus -> Aberdeen
     % H = interp1(x_SfA*1000,H_SfA,L,"linear","extrap");
    % dh_SfA = diff(H);
    
    dh = (H_2-H_1)/length(L);

    D_mm = 1000*D;% mm
    epsD = eps/D;

    % m_dot = rho_a*Q_a;

    P_a_kPa = P_a/1000;
    % Q_a_day = Q_a*24*3600;    % Sm3/day
    L_km = L./1000;           % km

    f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial friction f estimation using Nikuradse eq.

    rho = zeros(length(t),length(L));
    nu= zeros(length(t),length(L));
    nu_Po = zeros(length(t),length(L));
    Z = zeros(length(t),length(L));
    u = zeros(length(t),length(L));
    Re = zeros(length(t),length(L));
    f = zeros(length(t),length(L));
    s = zeros(length(t),length(L));
    h = zeros(length(t),length(L));
    psi = zeros(length(t),length(L));
    U_erosional = zeros(length(t),length(L));
    P = zeros(length(t),length(L));

    T = T_f*ones(length(t),length(L),1);      % Assuming isothermal flow
    
    P(:,1) = P_in; % Constant inlet pressure
    P(:,2) = P_a;

    for j = 1:length(t)
        
        for i = 1:length(L)
            rho(j,i) = CP.PropsSI('D','P',P(j,i),'T',T_f,'Air');
            nu(j,i) = CP.PropsSI('V','P',P(j,i),'T',T_f,'Air');
            h(j,i) = CP.PropsSI('H','P',P(j,i),'T',T_f,'Air');
            s(j,i) = CP.PropsSI('S','P',P(j,i),'T',T_f,'Air');
            Z(j,i) = CP.PropsSI('Z','P',P(j,i),'T',T_f,'Air');
            U_erosional(j,i) = 1.22*100/sqrt(rho(j,i)); % Erosional velocity - it is advised that the actual 
                                                    % velocity be up to 50% of the erosional velocity.
                                                    % The constant 100, can be any value from 100 to 250
    
            nu_Po(j,i) = 10*nu(j,i);    % Conversion of nu to Poise

            Q_a_day = 1.1494e-3*(T_a/P_a)*( (P(j,1)^2 - P(j,2)^2)/(G*T_f*dL*Z(j,i)*f_guess) )^0.5*D_mm^2.5;

            Q_old = Q_a_day;
            dQ = 10;
            count_Q = 0;
            
            while (dQ > 0.0001 & count_Q < 1000)
                Re(j,i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P in kPa, Q in m3/day, nu in poise (1 Pa s = 10 poise), and D in mm
            
                % Friction factor
                % Iterations for the estimation of friction factor using Colebrook equation
                f_old = f_guess;
                df = 10;
                count_f=0;
                while (df > 0.0001 & count_f < 10) 
                    % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
                    f_new = (-2*log10(epsD/3.7 + 2.825/(Re(j,i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
                    df = (f_new - f_old)/f_old;
                    f_old = f_new;
                    count_f = count_f + 1; 
                end
                f(j,i) = f_old;

                Q_a_day = 1.1494e-3*(T_a/P_a)*( (P(j,1)^2 - P(j,2)^2)/(G*T_f*dL*Z(j,i)*f(j,i)) )^0.5*D_mm^2.5;
                dQ = abs((Q_a_day - Q_old)/Q_old);
                Q_old = Q_a_day;
                count_Q = count_Q + 1;
            end
            
            Re(j,i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P in kPa, Q in m3/day, nu in poise (1 Pa s = 10 poise), and D in mm


            % Velocity m/s
            u(j,i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(j,i)*T_f/P(j,i)); % Q_b has to be converted to m3/day and D to mm
            if u(j,i) > 0.5*U_erosional(j,i)
                warning('U greater than 50% of the erosional velocity')
            end

            % Exergy
            % phi = h-h0 + T0*(s-s0)+u^2/2+gz
            % Considering height difference
            % psi(i) = h(i)-h0 + T0*(s(i)-s0)+u(i)^2/2+g*z(i)
            
            % Negligible height difference
            psi(j,i) = h(j,i)-h0 - T0*(s(j,i)-s0)+u(j,i)^2/2;
            
            







        end




        
    end

end
%%



for j=1:length(Var)
    % H_2 = Var{j};
    T_sin = Var{j};
    % dh_amp = Var{j};
    % H_prof = Var{j};
    % Flow_eq = Var{j};
    % D = Var{j};
    % Q_a = Var2{j};

    A = pi*D^2/4; % mm2
    L = 0:dL:L_m; % 70 km pipe with increments of 1 km

    % Elevation profile from St. Fergus -> Aberdeen
    % H = interp1(x_SfA*1000,H_SfA,L,"linear","extrap");
    % dh_SfA = diff(H);
    
    dh = (H_2-H_1)/length(L);

    D_mm = 1000*D;% mm
    epsD = eps/D;
    
    
    
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
    
    f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial friction f estimation using Nikuradse eq.

    for i=1:length(L)-1
        
        rho(i) = CP.PropsSI('D','P',P(i),'T',T_f,'Air');
        nu(i) = CP.PropsSI('V','P',P(i),'T',T_f,'Air');
        h(i) = CP.PropsSI('H','P',P(i),'T',T_f,'Air');
        s(i) = CP.PropsSI('S','P',P(i),'T',T_f,'Air');
        Z(i) = CP.PropsSI('Z','P',P(i),'T',T_f,'Air');
        U_erosional(i) = 1.22*100/sqrt(rho(i)); % Erosional velocity - it is advised that the actual 
                                                % velocity be up to 50% of the erosional velocity.
                                                % The constant 100, can be any value from 100 to 250

        nu_Po(i) = 10*nu(i);    % Conversion of nu to Poise
    




        % Velocity m/s
        u(i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(i)*T_f/P(i)); % Q_b has to be converted to m3/day and D to mm
        if u(i) > 0.5*U_erosional(i)
            warning('U greater than 50% of the erosional velocity')
        end

        % Exergy
        % phi = h-h0 + T0*(s-s0)+u^2/2+gz
        % Considering height difference
        % psi(i) = h(i)-h0 + T0*(s(i)-s0)+u(i)^2/2+g*z(i)

        % Negligible height difference
        psi(i) = h(i)-h0 - T0*(s(i)-s0)+u(i)^2/2;
        
        Re(i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P in kPa, Q in m3/day, nu in poise (1 Pa s = 10 poise), and D in mm
    
        % Friction factor
        % Iterations for the estimation of friction factor using Colebrook equation
        f_old = f_guess;
        df = 10;
        count_f=0;
        while (df > 0.0001 & count_f < 10) 
            f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
            f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
            df = abs((f_new - f_old)/f_old);
            f_old = f_new;
            count_f = count_f + 1; 
        end
        f(i) = f_old;
    
        P_1_kPa = P(i)/1000;
        P_2_kPa = 0.98*P_1_kPa; % Initial guess
        dP = 10;
        count_f = 0;
        
        while (dP > 0.0001 & count_f < 10)
            P_f_kPa = 2/3*( P_1_kPa + P_2_kPa-(P_1_kPa*P_2_kPa)/(P_1_kPa+P_2_kPa) );
            Z_f = CP.PropsSI('Z','P',P_f_kPa*1000,'T',T_f,'Air'); % compressibility for air
            % Compressibility formula for natural gas
            % Z_f = 1/(1+
            % ((P_f*0.000145038-14.73)*344400*10^(1.785*G)/(T_f*1.8)^3.825));
            
            if (H_prof == "Horizontal")

                if Flow_eq == "GFE"
                    % P from General Flow Equation
                    % Negligible height difference
                    P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(dL/1000)*Z_f*f(i) );
                elseif Flow_eq == "PanB"
                    % P from Panhandle B equation
                    % Negligible height difference
                    P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(dL/1000)*Z_f );
                    
                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                else
                    % P from Panhandle B equation
                    % Negligible height difference
                    P(i+1) = sqrt( P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(dL/1000)*Z_f );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                end

            elseif (H_prof == "Fixed tilt")
                
                s_H = 0.0684*G*dh/(T_f*Z_f);
                % L_e = dL*(exp(s_H)-1)/s_H;
                L_e = (dL/cos(atan(dh/dL)))*(exp(s_H)-1)/s_H;
                
                if Flow_eq == "GFE"
                    % P from General Flow Equation
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(L_e/1000)*Z_f*f(i) )/exp(s_H) );
                elseif Flow_eq == "PanB"
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                else
                    % P from Panhandle B equation
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                end

            elseif H_prof == "Sinusoidal"
                
                dh_sin = dh+dh_amp*sin( 2*pi/T_sin*L(i) );
                s_H = 0.0684*G*dh_sin/(T_f*Z_f); % Sinusoidal profile
                % L_e = dL*(exp(s_H)-1)/s_H;
                L_e = (dL/cos(atan(dh_sin/dL)))*(exp(s_H)-1)/s_H;
                
                if Flow_eq == "GFE"
                    % P from General Flow Equation
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(L_e/1000)*Z_f*f(i) )/exp(s_H) );
                elseif Flow_eq == "PanB"
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                else
                    % P from Panhandle B equation
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                end

            elseif H_prof == "Custom prof"

                dh = dh_SfA(i); % Custom elevation profile
        
                s_H = 0.0684*G*dh/(T_f*Z_f);
                % L_e = dL*(exp(s_H)-1)/s_H;
                L_e = (dL/cos(atan(dh/dL)))*(exp(s_H)-1)/s_H;

                if Flow_eq == "GFE"
                    % P from General Flow Equation
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.1494e-3*(T_a/P_a_kPa)*D_mm^2.5) )^2*G*T_f*(L_e/1000)*Z_f*f(i) )/exp(s_H) );
                elseif Flow_eq == "PanB"
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                else
                    % P from Panhandle B equation
                    P(i+1) = sqrt( (P_1_kPa^2 - ( Q_a_day/(1.002e-2*E*(T_a/P_a_kPa)^1.02*D_mm^2.53) )^(1/0.51)*G^0.961*T_f*(L_e/1000)*Z_f )/exp(s_H) );

                    % Panhandle B equation - valid for large diameter, high pressure flows with 4M < Re < 40M
                    if (Re(i) < 4e6) || (Re(i) > 40e6)
                        warning('Re outside the indicated region for the Panhandle B equation.')
                    end
                end
            else
                warning("Unidentified height profile.")

            end
            
            % CHECK CONDITIONS OF EQUATIONS

            dP = abs((P(i+1) - P_2_kPa)/P_2_kPa);
            P_2_kPa = P(i+1);
            count_f = count_f + 1;
        end
        P(i+1) = 1000*P(i+1); % conversion back to Pa
    
    end
    
    rho(end) = CP.PropsSI('D','P',P(end),'T',T_f,'Air');
    h(end) = CP.PropsSI('H','P',P(end),'T',T_f,'Air');
    s(end) = CP.PropsSI('S','P',P(end),'T',T_f,'Air');
    nu(end) = CP.PropsSI('V','P',P(end),'T',T_f,'Air');
    Z(end) = CP.PropsSI('Z','P',P(end),'T',T_f,'Air');
    u(end) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(end)*T_f/P(end)); % Q_b has to be converted to m3/day and D to mm
    psi(end) = h(end)-h0 - T0*(s(end)-s0)+u(end)^2/2;
    U_erosional(end) = 1.22*100/sqrt(rho(end)); % Erosional velocity - it is advised that the actual 
                                                % velocity be up to 50% of the erosional velocity.
                                                % The constant 100, can be any value from 100 to 250
    
    % Plots
    
    % Individual plots
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
    plot(L/1000,m_dot*psi/1e9,LStyle{j})
    xlabel('L [km]')
    ylabel('\Psi [GW]')
    % Percentage loss per km
    figure(pct_Psi_fig)
    Psi = m_dot*psi/1e9;
    Psi_n = 100*Psi./max(Psi);
    plot(L./1000,Psi_n,LStyle{j})    
    xlabel('L [km]')
    ylabel('\Psi/\Psi_{max} [%]')

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
    if H_prof == "Horizontal"
        plot(L/1000,zeros(length(L),1),LStyle{j})
    elseif H_prof == "Fixed tilt"
        plot(L/1000, (1:length(L))*dh,LStyle{j})
    elseif H_prof == "Sinusoidal"
        plot(L/1000, (0:length(L)-1)*dh + dh_amp*sin(2*pi/T_sin*L),LStyle{j})
    elseif H_prof == "Custom prof"
        plot(L/1000, H,LStyle{j})
    else
        plot(L/1000, (1:length(L))*dh,LStyle{j})
    end
    
    xlabel('L [km]')
    ylabel('H [m]')

end

figure(P_fig)
legend(string(Var))
applystyle2plot
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

figure(Psi_fig)
legend(string(Var))
applystyle2plot

% Percentage loss per km
figure(pct_Psi_fig)
legend(string(Var))
applystyle2plot
% Psi = m_dot*psi/1e9;
% Psi_n = Psi./max(Psi);
% plot(L./1000,Psi_n)
% grid on
% xlabel('L [km]')
% ylabel('\Psi/\Psi_{max}')

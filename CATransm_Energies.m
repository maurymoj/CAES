% pressure over 70 km pipeline
pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter

%% Steady-state analysis with constant T - X km pipe - Panhandle B flow equation
clear
clc
% close all
CP = py.importlib.import_module('CoolProp.CoolProp');
%--------------------------- Setup ---------------------------------



%-------------------- Problem parameters ------------------------------%
% Ambient conditions
T_a = 15 + 273.15; % T = 15 oC - ~288 K Standard temp used by the UK national gas company
P_a = 101325 ; % P = 101.325 kPa
G = 1;      % Specific gas gravity - for air G = 1

P_in = 9e6; % 7 MPa
% T_1 = 5 + 273.15; % Common operational condition assumption Nasr and Connor "Natural Gas Engineering and Safety Challenges"
T_1 = T_a; 
T_f = T_1; % Isothermal assumption
    % !!! High temperatures on the outlet of compressor stations,
    % which can persist for up to 50 km. !!!

% Flow rate
Q_a = 14*1000000/(24*3600); % Conversion from mcmd (Millions of standard cubic meters
    % per day [Mscmd] to standard cubic meters per second scms) - real demand estimation 70 mscm/day
    % St Fergus, 14 is the proportional relative to area of one of the 3, 900
    % mm diameter pipes
% Q_a = 25*1000000/(24*3600); % Flow for 1200 mm pipe
% Q_a = 3*1000000/(24*3600);  % Flow for 450 mm pipe

% Pipe
D = 0.9; % 900 mm
% D = 1.2;
% D = 0.45;
dL = 1000;    % Distance increment
L_m = 100000; % Total distance [m]
% L_m = 395000; % Total distance [m]
% L_m = 72145; % Total distance SfA [m]
% L_m = 73514; % Total distance SfA2 [m]


T_turb = 250+273.15;


eps = 0.04e-3; % Absolute roughness 0.04 mm

E = 0.75;   % Pipeline efficiency
    % Usually varies from 0.6 to 0.92 as a function of liquid presence and age
    % "Handbook of natural gas transmission and processing"

H_1 = 0; % Elevation at pipeline inlet
H_2 = 1; % Elevation at pipeline outlet

% H_1 = 30.6726; % Values from St. Fergus to Aberdeen profile 1
% H_2 = 137.917;

% H_1 = 8.1358; % Values from St. Fergus to Aberdeen profile 2
% H_2 = 104.6626;

T_sin = 25000; % Period of sinusoidal height increment
dh_amp = 25;    % Amplitude of sinusoidal height increment

% Elevation profile
% H_prof = "Horizontal";
H_prof = "Fixed tilt";
% H_prof = "Sinusoidal";
% H_prof = "Custom_prof";

% Flow equation

% Flow_eq = "GFE";
Flow_eq = "PanB"; % Overestimate of the pressure loss (~ 0.4% extra drop per 70 km)

% Generate multiple plots while varying one variable var (can be any
% parameter set in the initial statements)
% Pipe pressure Pa
% Var = {7e6, 8e6, 9e6, 10e6};
% Pipe diameter mm
% Var = {0.450,0.900,1.200}; Q = {3*1000000/(24*3600), 14*1000000/(24*3600), 25*1000000/(24*3600) };
% Var = {0.450,0.900,1.200}; Q = {1.9*3*1000000/(24*3600), 1.9*14*1000000/(24*3600), 1.9*25*1000000/(24*3600) };
% Var = {0.900,1.200};
% Var = {0.9}; Q = {15*1000000/(24*3600)};
% Flow rate scms (Standard cubic meters per second)
% Var2 = {3*1000000/(24*3600),14*1000000/(24*3600),25*1000000/(24*3600)};
% Var = {3*1000000/(24*3600),14*1000000/(24*3600),25*1000000/(24*3600)};
% Flow_rates = (3:25).*1000000/(24*3600);
% Elevation at H2
% Var = {100};
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


diams = [0.45, 0.9,1.2];

% Flow rate figure
% Flow_rates = (1:25)*1000000/(24*3600);

% Pressure figure
Ps = 1e6*(4:0.25:10);
Flow_rates = [3, 14, 25]*1000000/(24*3600); % 3, 14, and 25 Mscmd (Million of std m3 per day)

LStyle = {'b','r','k','b--','r--','k--','b-.'};

% % Sanity check with book equation comparisons
% G = 0.6
% T_a = (60 - 32)*5/9 + 273.15 % T = 15 oC - ~288 K
% P_a = 14.73*6894.76 % P = 14.73 psia
% rho_a = CP.PropsSI('D','P',P_a,'T',T_a,'Air');
% L_m = 161000                           % m
% D = (16 - 2*0.250)*0.0254              % m
% Q_a = (100*0.02825)*1000000/(24*3600)     % Sm3/s
% T_f = (80 - 32)*5/9 + 273.15
% P_in = (1400 + 14.73)*6894.76
% eps = 0.7e-3*25.4e-3
% E = 0.95

%----------------------------------------------------------------------%

rho_a = CP.PropsSI('D','P',P_a,'T',T_a,'Air');
T0 = T_a;
h0 = CP.PropsSI('H','P',P_a,'T',T_a,'Air');
s0 = CP.PropsSI('S','P',P_a,'T',T_a,'Air');

cp = CP.PropsSI('C','P',P_a,'T',T_a,'Air');
cv = CP.PropsSI('Cvmass','P',P_a,'T',T_a,'Air');
% gam = cp/cv;

% load('StF_Ab.mat');
load('StF_Ab2.mat');

Var1 = Ps;
Var2 = diams;

psi_dot_L = zeros(length(Var1),length(Var2));
dpsi_dot = zeros(length(Var1),length(Var2));

for k = 1:length(Var2)
    D = diams(k);
    Q_a = Flow_rates(k);
    
    for j = 1:length(Var1)
        % Inlet pressure
        % P_in = Var{j};
        % H_2 = Var{j};
        % T_sin = Var{j};
        % dh_amp = Var{j};
        % H_prof = Var{j};
        % Flow_eq = Var{j};    
        % Varying diameter and setting a proportional flow rate 
        % D = Var{j}; Q_a = Q{j}; 
        % Varying flow rate [scms]
        % Q_a = Flow_rates(j);
        P_in = Var1(j);
    
        A = pi*D^2/4; % mm2
        L = 0:dL:L_m; % 70 km pipe with increments of 1 km
    
        % Elevation profile from St. Fergus -> Aberdeen
        H = interp1(x_SfA*1000,H_SfA,L,"linear","extrap");
        dh_SfA = diff(H);
        
        dh = (H_2-H_1)/length(L);
    
        D_mm = 1000*D;% mm
        epsD = eps/D;
        
        m_dot = rho_a*Q_a;
    
    
        W_C = m_dot*cp*T_a*( (P_in/P_a)^(2/7) - 1);
    
    
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
    
        W_T = zeros(length(L),1);
        W_C2 = zeros(length(L),1);
        
        P(1) = P_in;
        
        f_guess = (2*log10(1/epsD)+1.14)^(-2); % Initial friction f estimation using Nikuradse eq.
    
        for i=1:length(L)-1
            
            rho(i) = CP.PropsSI('D','P',P(i),'T',T_f,'Air');
            nu(i) = CP.PropsSI('V','P',P(i),'T',T_f,'Air');
            h(i) = CP.PropsSI('H','P',P(i),'T',T_f,'Air');
            s(i) = CP.PropsSI('S','P',P(i),'T',T_f,'Air');
            Z(i) = CP.PropsSI('Z','P',P(i),'T',T_f,'Air');
    
            W_T(i) = m_dot*cp*T_turb*(1 - (P_a/P(i))^(2/7));
            W_C2(i) = m_dot*cp*T_f*((P_in/P(i))^(2/7) - 1);
    
            U_erosional(i) = 1.22*100/sqrt(rho(i)); % Erosional velocity - it is advised that the actual 
                                                    % velocity be up to 50% of the erosional velocity.
                                                    % The constant 100, can be any value from 100 to 250
    
            nu_Po(i) = 10*nu(i);    % Conversion of nu to Poise
        
            % Velocity m/s
            u(i) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(i)*T_f/P(i)); % Q_b has to be converted to m3/day and D to mm
            if u(i) > 0.5*U_erosional(i)
                break;
                warning('U greater than 50% of the erosional velocity')
            end
    
            % Exergy
            % phi = h-h0 + T0*(s-s0)+u^2/2+gz
    
            % Exergy considering height difference
            % psi(i) = h(i)-h0 + T0*(s(i)-s0)+u(i)^2/2+g*z(i)
            
            % Exergy assuming constant temperature
            % Negligible height difference
            psi(i) = h(i)-h0 - T0*(s(i)-s0)+u(i)^2/2;
            
            Re(i) = 0.5134*( P_a_kPa/T_a )*( G*Q_a_day/(nu_Po(i)*D_mm) ); % P in kPa, Q in m3/day, nu in poise (1 Pa s = 10 poise), and D in mm
        
            % Friction factor
            % Iterations for the estimation of friction factor using Colebrook equation
            f_old = f_guess;
            df = 10;
            n_iter=0;
            while (df > 0.0001 & n_iter < 10) 
                % f_n = (-2*log10(epsD/3.7 + 2.51/(Re*f^0.5)))^(-2); % Original Colebrook-White equation
                f_new = (-2*log10(epsD/3.7 + 2.825/(Re(i)*f_old^0.5)))^(-2); % Modified Colebrook-White equation - conservative (higher friction factors)
                df = abs((f_new - f_old)/f_old);
                f_old = f_new;
                n_iter = n_iter + 1; 
            end
            f(i) = f_old;
        
            P_1_kPa = P(i)/1000;
            P_2_kPa = 0.98*P_1_kPa; % Initial guess
            dP = 10;
            n_iter = 0;
            
            while (dP > 0.0001 & n_iter < 10)
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
                n_iter = n_iter + 1;
            end
            P(i+1) = 1000*P(i+1); % conversion back to Pa
            
            
            
    
        end
        
        if u(i) <= 0.5*U_erosional(i)
            rho(end) = CP.PropsSI('D','P',P(end),'T',T_f,'Air');
            h(end) = CP.PropsSI('H','P',P(end),'T',T_f,'Air');
            s(end) = CP.PropsSI('S','P',P(end),'T',T_f,'Air');
            nu(end) = CP.PropsSI('V','P',P(end),'T',T_f,'Air');
            Z(end) = CP.PropsSI('Z','P',P(end),'T',T_f,'Air');
        
            W_T(end) = m_dot*cp*T_turb*(1 - (P_a/P(end))^(2/7));
            W_C2(end) = m_dot*cp*T_f*((P_in/P(end))^(2/7) - 1);
        
            u(end) = 14.7349*( Q_a_day/D_mm^2 )*( P_a/T_a )*(Z(end)*T_f/P(end)); % Q_b has to be converted to m3/day and D to mm
            psi(end) = h(end)-h0 - T0*(s(end)-s0)+u(end)^2/2;
            U_erosional(end) = 1.22*100/sqrt(rho(end)); % Erosional velocity - it is advised that the actual 
                                                        % velocity be up to 50% of the erosional velocity.
                                                        % The constant 100, can be any value from 100 to 250

            psi_dot_L(j,k) = m_dot*psi(end); % flow exergy
            dpsi_dot(j,k) = m_dot*( psi(end) - psi(1) )/psi_dot_L(j,k); % exergy loss
        else
            psi_dot_L(j,k) = NaN; % flow exergy
            dpsi_dot(j,k) = NaN; % exergy loss
         end



    end

    
    
end


% Flow rate figure
% figure('Color',[1 1 1]);
% yyaxis left
% plot(Flow_rates*(24*3600)./1000000, psi_dot_L(:,1)./1e6)
% hold on 
% for k=2:length(diams)
%     plot(Flow_rates*(24*3600)./1000000, psi_dot_L(:,k)./1e6)
% end
% ylabel('$\dot \Psi$ [MW]','Interpreter','latex')
% % legend(num2str(diam'))
% 
% yyaxis right
% plot(Flow_rates*(24*3600)./1000000, dpsi_dot(:,1))
% hold on
% for k=2:length(diams)
%     plot(Flow_rates*(24*3600)./1000000, dpsi_dot(:,k))
% end
% ylabel('$\Delta \dot\Psi/\Psi_{in}$','Interpreter','latex')
% legend(num2str(diams'))
% xlabel('$Q_{st}$ [Mscmd]','Interpreter','latex')
% grid on

% Pressure figure
figure('Color',[1 1 1]);
yyaxis left
plot(Ps./1000000, psi_dot_L(:,1)./1e6)
hold on 
for k=2:length(diams)
    plot(Ps./1000000, psi_dot_L(:,k)./1e6)
end
ylabel('$\dot \Psi$ [MW]','Interpreter','latex')
% legend(num2str(diam'))

yyaxis right
plot(Ps./1000000, dpsi_dot(:,1))
hold on
for k=2:length(diams)
    plot(Ps./1000000, dpsi_dot(:,k))
end
ylabel('$\Delta \dot\Psi/\Psi_{in}$','Interpreter','latex')
legend(num2str(diams'))
xlabel('$P_{in}$ [MPa]','Interpreter','latex')
grid on
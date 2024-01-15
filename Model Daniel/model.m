%% St. Fergus to Aberdeen
%{
This simpler model will take into account the charging and discharging
capabilitites of a 2-stage compressor, 2 stage expansion ACAES installed at
St. Fergus gas terminal (scotland), and the pipeline connecting it to
Aberdeen, in Carlogie Gas compressor station. 
Info regarding the pipe network is taken from gas transmission network
files.

    -4 lines connect St. Fergus to Aberdeen directly. 

    Line        Diameter(m)         Length (m)      Approx. length (km)
     1            0.450              77 838               80
     2            0.900              64 089               65
     3            0.900              63 836               65
     4            0.900              67 156               65

Assumptions are:
(1) The pies are mainained at constant temperature (15C) by the ground acting a heat
sink/source whenever necessary
(2) Constant heat transfer coefficient between compressed gas and pipe,
30W/m²K
(3) each pipe will be discretised into 500m length sections

(4) The compressor and expander will operate with constant mass flow rate
and isentropic efficiency
(5) The pipe will be charged from 50 to 70 bar

%}


%% Loading CP if necessary
clearvars -except CP
if exist('CP','var')
else
    CP = py.importlib.import_module('CoolProp.CoolProp'); %If coolprop not loaded, then load it
end
wf = 'AIR';
%% Ambient conditions
T_a = 273.15 + 15; %K, Ambient temperature
P_a = 101325; %Pa, ambient pressure
u_a = CP.PropsSI('UMASS','P',P_a,'T',T_a,wf); %J/kg, air specific internal energy
h_a = CP.PropsSI('HMASS','P',P_a,'T',T_a,wf); %J/kg, air specific enthalpy
s_a = CP.PropsSI('SMASS','P',P_a,'T',T_a,wf); %J/kgK, air specific entropy
rho_a = CP.PropsSI('DMASS','P',P_a,'T',T_a,wf); %kg/m^3, air density
v_a = 1/rho_a; %m^3/kg, air specific volume

alpha_a = u_a + P_a*v_a - T_a*s_a; %specific non-flow exergy alpha function
beta_a = h_a - T_a*s_a; %air specific flow exergy beta function



%% Pipe conditions
d_pipe = [0.450;0.900;0.900;0.900]; %m, pipes diameters
L_total_pipe = [80;65;65;65]*1e3;%m, pipe approximated length
A_cs_pipe = pi.*d_pipe.^2/4; %m^2, cross section area
A_ht_pipe = pi.*d_pipe.*L_total_pipe; %m^2, total heat transfer area
V_pipe = pi.*d_pipe.^2/4.*L_total_pipe; %m^3, pipe volume
V_storage_total = sum(V_pipe);  %Total storage volume
T_pipe = T_a; %storage temperature
h_conv = 30; %W/m2K, heat transfer coefficient

%% Storage conditions

P_min = 50*1e5; %Pa, minimum storage pressure (50bar)
P_max = 70*1e5; %Pa, Maximum storage pressure (70bar)
P_thv = P_min; %Pa, throttling pressure

%% Heat exchangers
%Heat exchangers are assumed balanced Cr=1 for simplicity

delta_P_HX = 0.95; %5% pressure drop in heat exchangers
epsilon = 0.90;
T_cool_in = T_a;
cp_cool = 2000; %J/kgK Average coolant cp, loosely based on Dowtherm Q fluid


%% Flowmachinery

m_dot = 100; %kg/s, charging and discharging mass flow rate
eta_c = 0.85; %compressor isentropic efficiency
N_c = 2; %number of compressor units
eta_t = 0.90; %Expander isentropic efficiency
N_e = 2; %number of expander units
PR_t = (P_thv/P_a)^(1/N_e)*delta_P_HX;  %expander expansion ratio (constant due to throttle)

%% Cells
    %HPST storage
cell_storage = cell(2,6);
cell_storage(1,1:6) = {'SR','Pressure (Pa)','Temperature (K)','Mass (kg)','Non-flow Exergy (J)','Heat exchanged (W)'};
for m=2:6
    for n=2:6
        cell_storage{m,1} = m-1;
        cell_storage{m,n}(1,1:4) = {'Charging','Idle','Discharging','Recovery'}; %
    end
end
clear m n
    %Compressors
cell_compressors = cell(2,6);
cell_compressors(1,1:end) = {'SR','Power (W)','T_out (K)','P_out (Pa)','Irr (W)', 'Cumul. energy (J)'};
for m=2:6
    for n=2:6
        cell_compressors{m,1} = m-1;
        cell_compressors{m,n}(1,1:2) = {'C1','C2'};
    end
end
clear m n

    %Expanders
cell_expanders = cell(2,6);
cell_expanders(1,1:end) = {'SR','Power (W)','T_out (K)','P_out (Pa)','Irr (W)', 'Cumul. energy (J)'};
for m=2:6
    for n=2:6
        cell_expanders{m,1} = m-1;
        cell_expanders{m,n}(1,1:2) = {'E1','E2'};
    end
end
clear m n

    %Heat exchangers
cell_HX_chg = cell(2,7);
cell_HX_chg(1,1:end) = {'SR','Heat (W)','T_out air (K)','P_out air  (Pa)','T_cool_out  (K)','m_dot_cool (kg/s)','Irr  (W)'};

for m=2:6
    for n=2:7
        cell_HX_chg{m,1} = m-1;
        cell_HX_chg{m,n}(1,1:2) = {'IC','AC'};
    end
end
clear m n

cell_HX_dis = cell(2,7);
cell_HX_dis(1,1:end) = {'SR','Heat (W)','T_out air (K)','P_out air  (Pa)','T_cool_out  (K)','m_dot_cool (kg/s)','Irr  (W)'};

for m=2:6
    for n=2:7
        cell_HX_dis{m,1} = m-1;
        cell_HX_dis{m,n}(1,1:2) = {'HE1','HE2'};
    end
end
clear m n


    %TES
cell_TES_chg = cell(2,3);
cell_TES_chg(1,1:end) ={'SR','Mass (kg)','Temperature (K)'};
for m=2:6
    for n=2:3
        cell_TES_chg{m,1} = m-1;
        cell_TES_chg{m,n}(1,1:2) = {'ICTES','ACTES'};
    end
end
clear m n

cell_TES_dis = cell(2,3);
cell_TES_dis(1,1:end) ={'SR','Mass (kg)','Temperature (K)'};
for m=2:6
    for n=2:3
        cell_TES_dis{m,1} = m-1;
        cell_TES_dis{m,n}(1,1:2) = {'ICTES','ACTES'};
    end
end
clear m n

    % Throttle valve
cell_THV = cell(2,3);
cell_THV(1,1:end) = {'SR','T_out (K)','Irr (W)'};

for m=2:6
        cell_THV{m,1} = m-1;
end
clear m n
%% Iterative solution

SR=1; %initiating solution run

dt = 60; %1min timestep resolution
while SR<3  %Solution run indext, corresponds to the wait until isochoric system stabilises

clear time_*

    %% Charging
    clear t; t=1; %restarting time index
    time_chg(t,1) = 0;
    if SR==1 %first run - set values
%         P_st_chg(t,1) = P_min;
        cell_storage{SR+1,2}{2,1}(t,1) = P_min;
        cell_storage{SR+1,3}{2,1}(t,1) = T_a;
        cell_storage{SR+1,4}{2,1}(t,1) = CP.PropsSI('DMASS','P',cell_storage{SR+1,2}{2,1}(t,1),'T',cell_storage{SR+1,3}{2,1}(t,1),wf)*V_storage_total;
%         temp = (CP.PropsSI('UMASS','T',cell_storage{SR+1,3}{2,1}(t,1),'P',cell_storage{SR+1,2}{2,1}(t,1),wf) + P_a-CP.PropsSI('DMASS','T',cell_storage{SR+1,3}{2,1}(t,1),'P',cell_storage{SR+1,2}{2,1}(t,1),wf) - T_a*CP.PropsSI('SMASS','T',cell_storage{SR+1,3}{2,1}(t,1),'P',cell_storage{SR+1,2}{2,1}(t,1),wf));
        cell_storage{SR+1,5}{2,1}(t,1) = cell_storage{SR+1,4}{2,1}(t,1)*((CP.PropsSI('UMASS','T',cell_storage{SR+1,3}{2,1}(t,1),'P',cell_storage{SR+1,2}{2,1}(t,1),wf) + P_a-CP.PropsSI('DMASS','T',cell_storage{SR+1,3}{2,1}(t,1),'P',cell_storage{SR+1,2}{2,1}(t,1),wf) - T_a*CP.PropsSI('SMASS','T',cell_storage{SR+1,3}{2,1}(t,1),'P',cell_storage{SR+1,2}{2,1}(t,1),wf)) - alpha_a );
        
        
    else %not first run - load from previous SR
        cell_storage{SR+1,2}{2,1}(t,1) = cell_storage{SR,2}{2,3}(end,1); %Final recovery pressure on the previous solution run
        cell_storage{SR+1,3}{2,1}(t,1) = cell_storage{SR,3}{2,3}(end,1); %Final recovery temperature on the previous solution run
       cell_storage{SR+1,4}{2,1}(t,1) = cell_storage{SR,4}{2,3}(end,1); %Final recovery mass on the previous solution run
       cell_storage{SR+1,5}{2,1}(t,1) = cell_storage{SR,5}{2,3}(end,1);
    end

    %Resetting TES
        cell_TES_chg{SR+1,3}{2,1}(t,1) = 0;  %Initial temperature in ICTES = 0 
        cell_TES_chg{SR+1,3}{2,2}(t,1) = 0;  %Initial temperature in ACTES = 0 
        cell_TES_chg{SR+1,2}{2,1}(t,1) = 0;  %Initial mass stored in ICTES = 0 
        cell_TES_chg{SR+1,2}{2,2}(t,1) = 0;  %Initial mass stored in ACTES = 0 


    while cell_storage{SR+1,2}{2,1}(end,1)<P_max

    %Compressor 1
        PR_c(t,1) = (cell_storage{SR+1,2}{2,1}(t,1)/P_a)^(1/N_c)/delta_P_HX; %pressure ratio imposed on the compressor at time t
        [cell_compressors{SR+1,2}{2,1}(t,1),cell_compressors{SR+1,4}{2,1}(t,1),cell_compressors{SR+1,3}{2,1}(t,1),cell_compressors{SR+1,5}{2,1}(t,1)]...
            = compressor (m_dot, P_a, T_a, PR_c(t,1), eta_c, wf,CP,beta_a,T_a);

    % Intercooler
%         [Q_AHT,P_out_air,T_out_air,m_dot_cool,T_out_cool,Irr]...
          [cell_HX_chg{SR+1,2}{2,1}(t,1),cell_HX_chg{SR+1,4}{2,1}(t,1),cell_HX_chg{SR+1,3}{2,1}(t,1),cell_HX_chg{SR+1,6}{2,1}(t,1),cell_HX_chg{SR+1,5}{2,1}(t,1),cell_HX_chg{SR+1,7}{2,1}(t,1)]...
            = HX_chg(m_dot,cell_compressors{SR+1,3}{2,1}(t,1),cell_compressors{SR+1,4}{2,1}(t,1),cp_cool,T_a,epsilon,1,delta_P_HX,CP,beta_a,T_a,wf);

    % Compressor 2

 [cell_compressors{SR+1,2}{2,2}(t,1),cell_compressors{SR+1,4}{2,2}(t,1),cell_compressors{SR+1,3}{2,2}(t,1),cell_compressors{SR+1,5}{2,2}(t,1)]...
            = compressor (m_dot, cell_HX_chg{SR+1,4}{2,1}(t,1), cell_HX_chg{SR+1,3}{2,1}(t,1), PR_c(t,1), eta_c, wf,CP,beta_a,T_a);

    %Aftercooler
  %
       [cell_HX_chg{SR+1,2}{2,2}(t,1),cell_HX_chg{SR+1,4}{2,2}(t,1),cell_HX_chg{SR+1,3}{2,2}(t,1),cell_HX_chg{SR+1,6}{2,2}(t,1),cell_HX_chg{SR+1,5}{2,2}(t,1),cell_HX_chg{SR+1,7}{2,2}(t,1)]...
            = HX_chg(m_dot,cell_compressors{SR+1,3}{2,2}(t,1),cell_compressors{SR+1,4}{2,2}(t,1),cp_cool,T_a,epsilon,1,delta_P_HX,CP,beta_a,T_a,wf);
    

    %Storage update
% m_fin,T_fin,P_fin,Xi_fin = HPST(m_dot,T_flow,P_flow,m_ini,T_ini,P_ini,A_ht,V_st,h_conv,T_pipe,alpha_a,T_a,P_a,dt,CP,wf )
[cell_storage{SR+1,4}{2,1}(t+1,1),cell_storage{SR+1,3}{2,1}(t+1,1),cell_storage{SR+1,2}{2,1}(t+1,1),cell_storage{SR+1,5}{2,1}(t+1,1),cell_storage{SR+1,6}{2,1}(t+1,1)] ...
    = HPST(m_dot,cell_HX_chg{SR+1,3}{2,2}(t,1),cell_HX_chg{SR+1,4}{2,2}(t,1),cell_storage{SR+1,4}{2,1}(t,1),cell_storage{SR+1,3}{2,1}(t,1), cell_storage{SR+1,2}{2,1}(t,1),...
    sum(A_ht_pipe),V_storage_total,h_conv,T_pipe,alpha_a, T_a,P_a,dt,CP,wf);


    %TES update
[cell_TES_chg{SR+1,3}{2,1}(t+1,1),cell_TES_chg{SR+1,2}{2,1}(t+1,1),cell_TES_chg{SR+1,3}{2,2}(t+1,1),cell_TES_chg{SR+1,2}{2,2}(t+1,1)]...
    = TES_chg(cell_HX_chg{SR+1,5}{2,1}(t,1),cell_HX_chg{SR+1,6}{2,1}(t,1), cell_HX_chg{SR+1,5}{2,2}(t,1),cell_HX_chg{SR+1,6}{2,2}(t,1),cell_TES_chg{SR+1,3}{2,1}(t,1),cell_TES_chg{SR+1,2}{2,1}(t,1),cell_TES_chg{SR+1,3}{2,2}(t,1),cell_TES_chg{SR+1,2}{2,2}(t,1),cp_cool,dt);

    % Other updates
    time_chg(t+1,1) = time_chg(t,1)+dt; %updating time vector

    if t==1 %energy consumtion update
cell_compressors{SR+1,6}{2,1}(t,1) = 0 + cell_compressors{SR+1,2}{2,1}(t,1)*dt;  %Cumulative energy consumption at t=1 is 0 + power(t=1)*dt for C1
cell_compressors{SR+1,6}{2,2}(t,1) = 0 + cell_compressors{SR+1,2}{2,2}(t,1)*dt;  %Cumulative energy consumption at t=1 is 0 + power(t=1)*dt for C2
    else
cell_compressors{SR+1,6}{2,1}(t,1) = cell_compressors{SR+1,6}{2,1}(t-1,1) + cell_compressors{SR+1,2}{2,1}(t,1)*dt;  %Cumulative energy consumption at t=1 is cum. energy at t-1 + power(t=1)*dt for C1
cell_compressors{SR+1,6}{2,2}(t,1) = cell_compressors{SR+1,6}{2,2}(t-1,1) + cell_compressors{SR+1,2}{2,2}(t,1)*dt;  %Cumulative energy consumption at t=1 is cum. energy at t-1 + power(t=1)*dt for C2
    end

        t=t+1;

    end

time_chg=time_chg(1:end-1,1); %deleting last entry in time vector

    %% Idle
 %Idle will not be modelled now, as charging is basically isothermal.

%% Discharging

t=1; %resetting time index
time_dis(t,1)=0; %setting discharge time vector

    %Coupling discharge and charge processes
cell_storage{SR+1,2}{2,3}(t,1) = cell_storage{SR+1,2}{2,1}(end,1)  ; % Final charging pressure = initial discharging pressure
cell_storage{SR+1,3}{2,3}(t,1) = cell_storage{SR+1,3}{2,1}(end,1)  ; % Final charging temperature = initial discharging temperature
cell_storage{SR+1,4}{2,3}(t,1) = cell_storage{SR+1,4}{2,1}(end,1)  ; % Final charging stored mass = initial discharging stored mass
cell_storage{SR+1,5}{2,3}(t,1) = cell_storage{SR+1,5}{2,1}(end,1)  ; % Final charging stored non-flow exergy = initial discharging stored non-flow exergy 
cell_storage{SR+1,6}{2,3}(t,1) = cell_storage{SR+1,6}{2,1}(end,1)  ; % Final heat exchanged = initial heat exchanged
    
cell_TES_dis{SR+1,2}{2,1}(t,1) = cell_TES_chg{SR+1,2}{2,1}(end,1) ;  % Initial mass in ICTES discharging = final ICTES charging mass
cell_TES_dis{SR+1,3}{2,1}(t,1) = cell_TES_chg{SR+1,3}{2,1}(end,1) ;  % Initial temperature in ICTES discharging = final ICTES charging temperature 
cell_TES_dis{SR+1,2}{2,2}(t,1) = cell_TES_chg{SR+1,2}{2,2}(end,1) ;  % Initial mass in ACTES discharging = final ACTES charging mass
cell_TES_dis{SR+1,3}{2,2}(t,1) = cell_TES_chg{SR+1,3}{2,2}(end,1) ;  % Initial temperature in ACTES discharging = final ACTES charging temperature 


%iterative process
while cell_storage{SR+1,2}{2,3}(end,1)>P_min


    % Storage

[cell_storage{SR+1,4}{2,3}(t+1,1),cell_storage{SR+1,3}{2,3}(t+1,1),cell_storage{SR+1,2}{2,3}(t+1,1),cell_storage{SR+1,5}{2,3}(t+1,1),cell_storage{SR+1,6}{2,3}(t+1,1)] ...
    = HPST(-m_dot,cell_storage{SR+1,3}{2,3}(t,1),cell_storage{SR+1,2}{2,3}(t,1),cell_storage{SR+1,4}{2,3}(t,1),cell_storage{SR+1,3}{2,3}(t,1), cell_storage{SR+1,2}{2,3}(t,1),...
    sum(A_ht_pipe),V_storage_total,h_conv,T_pipe,alpha_a, T_a,P_a,dt,CP,wf);


    % Throttle valve

[cell_THV{SR+1,2}(t,1),cell_THV{SR+1,3}(t,1)] = throttle_valve(cell_storage{SR+1,2}{2,3}(t,1),cell_storage{SR+1,3}{2,3}(t,1),P_thv,m_dot,beta_a,T_a,wf,CP);
    
    % Heater 1
%[Q_AHT,P_out_air,T_out_air,m_dot_cool,T_out_cool,Irr] = HX_dis(m_dot_air,T_in_air,P_in_air,cp_cool,T_in_cool,epsilon,C_r,DP_air,CP,beta_a,T_a,wf)
[cell_HX_dis{SR+1,2}{2,1}(t,1),cell_HX_dis{SR+1,4}{2,1}(t,1),cell_HX_dis{SR+1,3}{2,1}(t,1),cell_HX_dis{SR+1,6}{2,1}(t,1),cell_HX_dis{SR+1,5}{2,1}(t,1),cell_HX_dis{SR+1,7}{2,1}(t,1)] = ...
    HX_dis(m_dot,cell_THV{SR+1,2}(t,1),P_thv,cp_cool,cell_TES_dis{SR+1,3}{2,1},epsilon,1,delta_P_HX,CP,beta_a,T_a,wf);


    % Expander 1
%[W_dot, P_out,T_out,Irr] = expander(m_dot,P_in,T_in,PR,eta_iso,wf,CP,beta_a,T_a)
    [cell_expanders{SR+1,2}{2,1}(t,1), cell_expanders{SR+1,4}{2,1}(t,1), cell_expanders{SR+1,3}{2,1}(t,1),cell_expanders{SR+1,5}{2,1}(t,1)] ...
        = expander(m_dot,cell_HX_dis{SR+1,4}{2,1}(t,1),cell_HX_dis{SR+1,3}{2,1}(t,1),1/PR_t,eta_t,wf,CP,beta_a,T_a);


 % Heater 2
%[Q_AHT,P_out_air,T_out_air,m_dot_cool,T_out_cool,Irr] = HX_dis(m_dot_air,T_in_air,P_in_air,cp_cool,T_in_cool,epsilon,C_r,DP_air,CP,beta_a,T_a,wf)
[cell_HX_dis{SR+1,2}{2,2}(t,1),cell_HX_dis{SR+1,4}{2,2}(t,1),cell_HX_dis{SR+1,3}{2,2}(t,1),cell_HX_dis{SR+1,6}{2,2}(t,1),cell_HX_dis{SR+1,5}{2,2}(t,1),cell_HX_dis{SR+1,7}{2,2}(t,1)] = ...
    HX_dis(m_dot,cell_expanders{SR+1,3}{2,1}(t,1),cell_expanders{SR+1,4}{2,1}(t,1),cp_cool,cell_TES_dis{SR+1,3}{2,2},epsilon,1,delta_P_HX,CP,beta_a,T_a,wf);


    % Expander 2
%[W_dot, P_out,T_out,Irr] = expander(m_dot,P_in,T_in,PR,eta_iso,wf,CP,beta_a,T_a)
    [cell_expanders{SR+1,2}{2,2}(t,1), cell_expanders{SR+1,4}{2,2}(t,1), cell_expanders{SR+1,3}{2,2}(t,1),cell_expanders{SR+1,5}{2,2}(t,1)] ...
        = expander(m_dot,cell_HX_dis{SR+1,4}{2,2}(t,1),cell_HX_dis{SR+1,3}{2,2}(t,1),1/PR_t,eta_t,wf,CP,beta_a,T_a);

 if t==1 %energy consumtion update
cell_expanders{SR+1,6}{2,1}(t,1) = 0 + cell_expanders{SR+1,2}{2,1}(t,1)*dt;  %Cumulative energy generation at t=1 is 0 + power(t=1)*dt for C1
cell_expanders{SR+1,6}{2,2}(t,1) = 0 + cell_expanders{SR+1,2}{2,2}(t,1)*dt;  %Cumulative energy generation at t=1 is 0 + power(t=1)*dt for C2
    else
cell_expanders{SR+1,6}{2,1}(t,1) = cell_expanders{SR+1,6}{2,1}(t-1,1) + cell_expanders{SR+1,2}{2,1}(t,1)*dt;  %Cumulative energy generation at t=1 is cum. energy at t-1 + power(t=1)*dt for C1
cell_expanders{SR+1,6}{2,2}(t,1) = cell_expanders{SR+1,6}{2,2}(t-1,1) + cell_expanders{SR+1,2}{2,2}(t,1)*dt;  %Cumulative energy generation at t=1 is cum. energy at t-1 + power(t=1)*dt for C2
    end
% Other updates
    time_dis(t+1,1) = time_dis(t,1)+dt; %updating time vector
t=t+1;
end
rte(SR,1) = (cell_expanders{2,6}{2,1}(end,1)+cell_expanders{2,6}{2,2}(end,1))/(cell_compressors{2,6}{2,1}(end,1)+cell_compressors{2,6}{2,2}(end,1));
SR=SR+1;
 end
% 

function [Q_AHT,P_out_air,T_out_air,m_dot_cool,T_out_cool,Irr] = HX_dis(m_dot_air,T_in_air,P_in_air,cp_cool,T_in_cool,epsilon,C_r,DP_air,CP,beta_a,T_a,wf)
%Basic heat exchanger calculation


%Inlet air conditions
h_in_air =  CP.PropsSI('HMASS','T', T_in_air, 'P',P_in_air,wf);  %Inlet specific enthalpy
cp_in_air = CP.PropsSI('CPMASS','T', T_in_air, 'P',P_in_air,wf);  %inlet specific cp
s_in_air =  CP.PropsSI('SMASS','T', T_in_air, 'P',P_in_air,wf);  %Inlet specific entropy

xi_in_air = h_in_air - T_a*s_in_air - beta_a;  % Inlet specific flow exergy
%Outlet air pressure
P_out_air = P_in_air*DP_air; %outlet pressure
%A heat transfer

h_MHT_out_air = CP.PropsSI('HMASS','P',P_out_air,'T',T_in_cool,wf); %outlet enthalpy at MHT conditions 
Q_AHT = epsilon*m_dot_air*(h_MHT_out_air-h_in_air); % Actual heat transferred

%Air outlet 
h_our_air = h_in_air + Q_AHT/m_dot_air; %outlet air specific enthalpy 
T_out_air = CP.PropsSI('T','P',P_out_air,'HMASS',h_our_air,wf);  %Outlet temperature
cp_air_out = CP.PropsSI('CPMASS','P',P_out_air,'HMASS',h_our_air,wf);  %Outlet cp
s_out_air = CP.PropsSI('SMASS','P',P_out_air,'HMASS',h_our_air,wf);  %Outlet specific entropy

%Coolant outlet

cp_air_mean = 0.5*(cp_in_air +cp_air_out); %
m_dot_cool = C_r*(m_dot_air*cp_air_mean)/(cp_cool); %From definition of heat capacity
T_out_cool = T_in_cool - Q_AHT/(m_dot_cool*cp_cool);


%% Exergy analysis

Irr = T_a*(m_dot_air*cp_air_mean*log(T_out_air/T_in_cool) - m_dot_cool*cp_cool*log(T_in_cool/T_out_cool) ) + m_dot_air*287*T_a*log(P_in_air/P_out_air) + m_dot_cool*cp_cool*(T_out_cool-T_a) ; 

end
function [T_out,Irr] = throttle_valve(P_in,T_in,P_out,m_dot,beta_a,T_a,wf,CP)
%Throttle valve
%The throttle valve function calculates the inlet enthalpy and exergy. The
%vale operation is isenthalpic, therefore, when expansion ration
%(P_thr/P_in<<1), a cooling effect is expected. the magnitude decreases as
%the discharging process take place.

%Inlet side

h_in = CP.PropsSI('HMASS','T',T_in,'P',P_in,wf); %inlet specific enthalpy
s_in = CP.PropsSI('SMASS','T',T_in,'P',P_in,wf); %inlet specific entropy

xi_in = (h_in - T_a*s_in)-beta_a; %inlet exergy flow

%Outlet side
h_out = h_in; %isenthalpic
T_out = CP.PropsSI('T','HMASS',h_out,'P',P_out,wf); %outlet temperature
s_out = CP.PropsSI('SMASS','HMASS',h_out,'P',P_out,wf); %outlet specific entropy

xi_out = (h_out - T_a*s_out)-beta_a; %inlet exergy flow

Irr = m_dot*abs(xi_out-xi_in); %irreversibility (effectively Guoy-Stodola since h_in = h_out) 
end 
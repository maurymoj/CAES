function [m_fin,T_fin,P_fin,Xi_fin,Q_st] = HPST(m_dot,T_flow,P_flow,m_ini,T_ini,P_ini,A_ht,V_st,h_conv,T_pipe,alpha_a,T_a,P_a,dt,CP,wf )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Calculating parameters

h_flow = CP.PropsSI('HMASS','T',T_flow,'P',P_flow,wf); % inflow (chg) or outflow (dis) air specific enthalpy
Q_st = -h_conv*A_ht*(T_ini-T_pipe); %heat exchanged

u_ini = CP.PropsSI('UMASS','T',T_ini,'P',P_ini,wf); %initial specific internal energy

%Balance of mass
m_fin = m_ini + m_dot*dt;

%Balance of energy
u_fin = ((Q_st+m_dot*h_flow)*dt+m_ini*u_ini)/(m_fin); %final specific internal energy


%Calculating updated properties
rho_fin = m_fin/V_st;  %updated stored air density 
T_fin = CP.PropsSI('T','DMASS',rho_fin,'UMASS',u_fin,wf); %updated temperature
P_fin = CP.PropsSI('P','DMASS',rho_fin,'UMASS',u_fin,wf); %updated pressure
s_fin = CP.PropsSI('SMASS','DMASS',rho_fin,'UMASS',u_fin,wf); %updated specific entropy
% Balance of exergy

Xi_fin = m_fin*((u_fin + P_a/rho_fin - T_a*s_fin)-alpha_a);



end
function [T_ICTES_fin, m_ICTES_fin,T_ACTES_fin, m_ACTES_fin ] = TES_chg(T_out_IC,m_dot_out_IC,T_out_AC,m_dot_out_AC,T_ICTES_ini,m_ICTES_ini,T_ACTES_ini,m_ACTES_ini,cp_cool,dt)
%TES_charging
% This function aaplies mass and energy conservation to the IC and AC TES,
% assuming constant cp_cool

%% For ICTES
%B.o.M
m_ICTES_fin = m_ICTES_ini + m_dot_out_IC*dt; %final mass in ICTES, after time increment
%B.o.E
T_ICTES_fin = (m_ICTES_ini*T_ICTES_ini*cp_cool + m_dot_out_IC*T_out_IC*cp_cool*dt)/(m_ICTES_fin*cp_cool); %Final mixed temperature

%% FOR ACTES
%B.o.M
m_ACTES_fin = m_ACTES_ini + m_dot_out_AC*dt; %final mass in ACTES, after time increment
%B.o.E
T_ACTES_fin = (m_ACTES_ini*T_ACTES_ini*cp_cool + m_dot_out_AC*T_out_AC*cp_cool*dt)/(m_ACTES_fin*cp_cool); %Final mixed temperature

end
function [W_dot, P_out,T_out,Irr] = expander(m_dot,P_in,T_in,PR,eta_iso,wf,CP,beta_a,T_a)
%expander function
%{
Inputs: m_dot,P_in,T_in,PR,eta_iso,wf,CP
Outputs: W_dot, P_out,T_out,Irr

Using definition of isentropic efficiency for expander, outlet
conditions are calculated.


%}
% Inlet conditions
h_in = CP.PropsSI('HMASS','T',T_in,'P',P_in,wf); % Inlet specific enthalpy
s_in = CP.PropsSI('SMASS','T',T_in,'P',P_in,wf); % Inlet specific entropy
xi_in = h_in - T_a*s_in - beta_a; %inlet specific flow exergy


%Outlet pressure
P_out = P_in*PR; %Pa, outlet pressure

%Isentropic conditions at outlet
h_outs = CP.PropsSI('HMASS','SMASS',s_in,'P',P_out,wf); % Outlet isentropic specific enthalpy
%Actual outlet
h_out = (h_outs-h_in)*eta_iso+h_in; %Outlet specific enthalpy
T_out = CP.PropsSI('T','HMASS',h_out,'P',P_out,wf); %outlet temperature
W_dot = m_dot*abs(h_out-h_in); % Power consumption
s_out = CP.PropsSI('SMASS','HMASS',h_out,'P',P_out,wf); %outlet specific entropy

xi_out = h_out - T_a*s_out - beta_a; %outlet specific flow exergy

% Exergy Balance

Irr = m_dot*T_a*(s_out-s_in); % Irreversibility



end
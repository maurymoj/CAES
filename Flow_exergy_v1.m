% Flow exergy

% Mass flow rate equation based on st. fergus data (5 pipelines, 1x d=0.4, 3x d=0.9 and 1x d=1.2m;
Total_area = pi*0.45^2/4 + 3*pi*0.9^2/4 + pi*1.2^2/4;
Total_flow = 42.6 + 3*198.6 + 354.6;


m_dot = 198.6; % [kg/s] - Value proportional to the area based on total mass flow rate into St. Fergus
P = 7e6; % 7 MPa
T = 298;

D = 0.900;

P_o = P_amb;
T_o = T_amb;
h_o = CP.PropsSI('H','P',P_o,'T',T_o,'air');
s_o = CP.PropsSI('S','P',P_o,'T',T_o,'air');

h = CP.PropsSI('H','P',P,'T',T,'air');
s = CP.PropsSI('S','P',P,'T',T,'air');
rho = CP.PropsSI('D','P',P,'T',T,'air');

A_h = pi*D^2/4;
V = m_dot/(rho*A_h);
psi =  (h - h_o) - T_o*(s - s_o)+V^2/2; % + g*z ;
Psi = m_dot*psi % [W]
Psi_MW = Psi/(1e6) % [MW]

CP = py.importlib.import_module('CoolProp.CoolProp'); % Simplifies coolprop calls

P = 0.1e6:0.1e6:7e6;
T = 278:328;

P_o = 101325;
T_o = 298.15;
rho_o = CP.PropsSI('D','P',P_o,'T',T_o,'air');
v_o = 1/rho_o;
u_o   = CP.PropsSI('U','P',P_o,'T',T_o,'air');
s_o   = CP.PropsSI('S','P',P_o,'T',T_o,'air');

rho = zeros(length(P),length(T));
u   = zeros(length(P),length(T));
s   = zeros(length(P),length(T));
v   = zeros(length(P),length(T));
chi = zeros(length(P),length(T));

for i=1:length(P)
    for j=1:length(T)
        rho(i,j) = CP.PropsSI('D','P',P(i),'T',T(j),'air');
        v(i,j) = 1/rho(i,j);
        u(i,j)   = CP.PropsSI('U','P',P(i),'T',T(j),'air');
        s(i,j)   = CP.PropsSI('S','P',P(i),'T',T(j),'air');
        chi(i,j) = ( (u(i,j) - u_o) + P_o*(v(i,j) - v_o) - T_o*(s(i,j) - s_o) );
    end
end

[X,Y] = meshgrid(P,T);
Z = chi;
figure('Color',[1 1 1])
surf(X,Y,Z')

CP = py.importlib.import_module('CoolProp.CoolProp');
rho_o = CP.PropsSI('D','P',P_o,'T',T_o,'air');
figure;plot(t./3600,m.*(mean(u)-u_o)/(3600*1e6),t./3600,P_o*m.*(1./mean(rho)-1/rho_o)/(3600*1e6),t./3600, -T_o*m.*(mean(s)-s_o)/(3600*1e6) )
title(heat_transfer_model)
grid on
ylabel('\chi [MWh]')

figure;plot(t./3600,X_coolprop - X_coolprop(1))
hold on
plot(t./3600,X - X(1))
grid on
ylabel('\Delta\chi [MWh]')
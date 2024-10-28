XX = sum(m_n(2:end,:).*( u(2:end,:) - u_o + P_o*R*(T(2:end,:)./P(2:end,:) - T_o/P_o) - T_o*(s(2:end,:) - s_o) ) )./(1e6*3600);

dXX = rho_f(2,:)'.*v(2,:)'*A_h.*(h_f(2,:)' - h_o - T_o*(s_f(2,:)' - s_o))*dt;

XX_bal = XX(1) + cumsum(dXX);

XX(end);
sum(dXX)/(1e6*3600);

figure
plot(t,XX,t,XX(1)+cumsum(dXX)./(1e6*3600))
legend('XX','XX_{bal}')

figure
plot((XX-XX_bal')./XX)
%%  Cycle analysis
figure
plot(t,XX,t,[XX(1)+cumsum(dXX(t<Dt_charg))./(1e6*3600);XX(j_disch)+cumsum(dXX(t>=Dt_charg))./(1e6*3600)])

figure
subplot(1,2,1)
plot(t(t<Dt_charg),XX(t<Dt_charg),t(t<Dt_charg),XX(1)+cumsum(dXX(t<Dt_charg))./(1e6*3600))
subplot(1,2,2)
plot(t(t>=Dt_charg),XX(t>=Dt_charg),t(t>=Dt_charg),XX(j_disch)+cumsum(dXX(t>=Dt_charg))./(1e6*3600))

figure('color',[1 1 1]);plot(t,m)
subplot(1,2,1)
t_charg = t(t<Dt_charg);
plot(t_charg,XX(t<Dt_charg),t_charg,XX(1)+cumsum(dXX(t<Dt_charg))./(1e6*3600))
title('Charging')
% ylim([2.9 5])
subplot(1,2,2)
t_disch = t(t>=Dt_charg) - t_charg(end);
plot(t_disch,XX(t>=Dt_charg),t_disch,XX(j_disch)+cumsum(dXX(t>=Dt_charg))./(1e6*3600))
title('Discharging')
% ylim([2.9 5])

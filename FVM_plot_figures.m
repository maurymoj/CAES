% Plot the full cycle of energy storage dividing the result in charging and
% discharging processes
if strcmp(Process,'Cycle_L') | strcmp(Process,'Cycle_R')
    figure('color',[1 1 1]);plot(t,m)
    hold on; plot(t,m_bal)
    legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
    figure('color',[1 1 1]);plot(t,E./(1e6*3600))
    hold on; plot(t,E_bal(1:end)./(1e6*3600))
    grid on 
    legend('E [MWh]','$E_o + \dot{m} \Delta E [MWh]$','Interpreter','latex')

    figure('color',[1 1 1]);
    subplot(1,2,1)
    t_charg = t(t<Dt_charg);
    plot(t_charg,X(t<Dt_charg),t_charg,X(1)+cumsum(dX(t<Dt_charg))./(1e6*3600))
    ylabel('Exergy')
    title('Charging')
    ylim([2.9 5])
    subplot(1,2,2)
    t_disch = t(t>=Dt_charg) - t_charg(end);
    plot(t_disch,X(t>=Dt_charg),t_disch,X(j_disch)+cumsum(dX(t>=Dt_charg))./(1e6*3600))
    ylabel('Exergy')
    title('Discharging')
    ylim([2.9 5])

    % figure('color',[1 1 1]);
    % subplot(1,2,1)
    % title('Charging')
    % plot(t(t<Dt_charg),X(t<Dt_charg),t(t<Dt_charg),X(1)+cumsum(dX(t<Dt_charg))./(1e6*3600))
    % subplot(1,2,2)
    % title('Discharging')
    % plot(t(t>=Dt_charg),X(t>=Dt_charg),t(t>=Dt_charg),X(j_disch)+cumsum(dX(t>=Dt_charg))./(1e6*3600))
    
else
    figure('color',[1 1 1]);plot(t,m)
    hold on; plot(t,m_bal)
    legend('m','$m_o + \dot{m} dt$','Interpreter','latex')
    figure('color',[1 1 1]);plot(t,E./(1e6*3600))
    hold on; plot(t,E_bal(1:end)./(1e6*3600))
    grid on
    legend('E [MWh]','$E_o + \dot{m} \Delta E [MWh]$','Interpreter','latex')
end

figure('Color',[1 1 1]);
plot(abs(mean(error_hist)))
xlabel('Iteration')
ylabel('Absolute mean of the error')
grid on

figure;
yyaxis left
plot(t,m2,t,m2_bal)
ylabel('mass [kg]')
yyaxis right
plot(t,E2./(1e6*3600),t,E2_bal./(1e6*3600))
ylabel('Energy [MWh]')
grid on
legend('m','m_{bal}','E','E_{bal}')

figure
plot((E2 - E2_bal')./E2)
xlabel(' t [s]')
ylabel('Relative difference between E and E_{bal}')
grid on

if strcmp(Process,'Kiuchi')
    figure('Color',[1 1 1])
    plot(t./60,rho_f(1,:).*v(1,:)*A_h*3600/rho_st)
    hold  on
    plot(t./60,rho_f(end,:).*v(end,:)*A_h*3600/rho_st)
    legend('Face 1','Face n')
    ylim([-1e5 4e5])
    ylabel('Volumetric flow rate [scmh]')

    figure('Color',[1 1 1])
    plot(t./60,rho(1,:).*v_n(1,:)*A_h*3600/rho_st,'b')
    hold  on
    plot(t./60,rho(end,:).*v_n(end,:)*A_h*3600/rho_st,'r')
    legend('Node 1','Node n')
    ylabel('Volumetric flow rate [scmh]')

    figure('Color',[1 1 1])
    plot(t./60,rho_f(1,:).*v(1,:)*A_h)
    hold  on
    plot(t./60,rho_f(end,:).*v(end,:)*A_h)
    legend('Face 1','Face n')
    ylabel('Mass flow rate [kg/s]')
end
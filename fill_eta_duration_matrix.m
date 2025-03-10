% Generates matrices with the values of storage efficiency, process
% duration (charging/discharging time), and elapsed time (total time for
% the conclusion of the simulation).

clear

Ps = 6:10;
Pmins = Ps - 3;
m_dots = 50:25:200;

% Process types:
% Charging_L
% Discharging_L
Process_type = 'Discharging_L';
etaX_stor_M = zeros(length(Ps),length(m_dots));
Duration = zeros(length(Ps),length(m_dots));
elapsedTimes = zeros(length(Ps),length(m_dots));

for ii=1:length(Ps)
    for jj = 1:length(m_dots)
        % filename = strcat('CAESPipe_Charging_L_P_',num2str(Ps(ii)),'MPa_m_',num2str(ms(jj)));
        % filename = strcat('CAESPipe_Charging_L_P_',num2str(Ps(ii)),'-',num2str(Pmins(ii)),'MPa_m_',num2str(ms(jj)));
        filename = strcat('CAESPipe_',Process_type,'_P_',num2str(Ps(ii)),'-',num2str(Pmins(ii)),'MPa_m_',num2str(m_dots(jj)))
        load(filename)
        % Ps(ii)
        % ms(jj)
        % t_shut_off
        etaX_stor_M(ii,jj) = etaX_stor;
        Duration(ii,jj) = t_shut_off;
        elapsedTimes(ii,jj) = elapsedTime;
        % figure('Color',[1 1 1]);
        % plot(abs(mean(error_hist)))
    end
end
% 
% figure('Color',[1 1 1])
% hold all
% for ii = 1:length(Ps)
%     plot(ms,elapsedTimes(ii,:))
% end
% grid on
% legend(num2str(Ps'))

%%
figure('Color',[1 1 1])
surf(m_dots,Ps,etaX_stor_M)
xlabel('$\dot{m}$ [kg/s]','Interpreter','latex')
ylabel('$P_{max}$ [MPa]','Interpreter','latex')
zlabel('$\eta_{Ch}$','Interpreter','latex')
title(strcat(Process,' efficiency [-]'))

figure('Color',[1 1 1])
hold all
for ii = 1:length(Ps)
    plot(m_dots,etaX_stor_M(ii,:))
end
grid on
xlabel('$\dot{m}$ [kg/s]','Interpreter','latex')
ylabel('$\eta_{Ch}$','Interpreter','latex')
legend(num2str(Ps'))
title(strcat(Process,' efficiency [-]'))


figure('Color',[1 1 1])
surf(m_dots,Ps,Duration./3600)
xlabel('m_{dot}')
ylabel('P_{max}')
zlabel(strcat(Process,' Duration [h]'))

figure('Color',[1 1 1])
hold all
for ii = 1:length(Ps)
    plot(m_dots,Duration(ii,:)./3600)
end
grid on
xlabel('$\dot{m}$ [kg/s]','Interpreter','latex')
ylabel(strcat(Process,' Duration [h]'))
legend(num2str(Ps'))
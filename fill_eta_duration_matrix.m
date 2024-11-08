clear

Ps = 6:10;
Pmins = Ps - 3;
ms = 50:50:200;
etaX_stor_M = zeros(length(Ps),length(ms));
Duration = zeros(length(Ps),length(ms));
elapsedTimes = zeros(length(Ps),length(ms));

for ii=1:length(Ps)
    for jj = 1:length(ms)
        % filename = strcat('CAESPipe_Charging_L_P_',num2str(Ps(ii)),'MPa_m_',num2str(ms(jj)));
        filename = strcat('CAESPipe_Charging_L_P_',num2str(Ps(ii)),'-',num2str(Pmins(ii)),'MPa_m_',num2str(ms(jj)));
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
surf(ms,Ps,etaX_stor_M)
xlabel('$\dot{m}$ [kg/s]','Interpreter','latex')
ylabel('$P_{max}$ [MPa]','Interpreter','latex')
zlabel('$\eta_{Ch}$','Interpreter','latex')

figure('Color',[1 1 1])
hold all
for ii = 1:length(Ps)
    plot(ms,etaX_stor_M(ii,:))
end
grid on
xlabel('$\dot{m}$ [kg/s]','Interpreter','latex')
ylabel('$\eta_{Ch}$','Interpreter','latex')
legend(num2str(Ps'))


figure('Color',[1 1 1])
surf(ms,Ps,Duration./3600)
xlabel('m_{dot}')
ylabel('P_{max}')
zlabel('Charging Duration [h]')

figure('Color',[1 1 1])
hold all
for ii = 1:length(Ps)
    plot(ms,Duration(ii,:)./3600)
end
grid on
xlabel('m_{dot} [kg/s]')
ylabel('Charging Duration [h]')
legend(num2str(Ps'))
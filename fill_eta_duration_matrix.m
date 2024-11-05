clear

Ps = 6:9;
Pmins = Ps - 3;
ms = 50:50:200;
% ms = 200;
etaX_stor_M = zeros(length(Ps),length(ms));
Duration = zeros(length(Ps),length(ms));

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
    end
end

figure('Color',[1 1 1])
surf(ms,Ps,etaX_stor_M)
xlabel('m_{dot} [kg/s]')
ylabel('P_{max} [MPa]')
zlabel('\eta_{ex}')

figure('Color',[1 1 1])
hold all
for ii = 1:length(Ps)
    plot(ms,etaX_stor_M(ii,:))
end
grid on
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
% Cavern x Pipeline
% Pressure profiles
% Charging

figure('color',[1 1 1])
colororder('sail')
load('CAESPipe_Charging_L_P_7-4MPa_m_100.mat')
plot(x_n./L,P(:,(t>=t_shut_off) & (t<t_shut_off+ dt))./1e6)
hold on
plot(x_n./L,P(:,(t>=t_shut_off + 0.5*3600) & (t<t_shut_off+ 0.5*3600+dt))./1e6)
plot(x_n./L,P(:,(t>=t_shut_off + 1*3600) & (t<t_shut_off+ 1*3600+dt))./1e6)
plot(x_n./L,P(:,(t>=t_shut_off + 1.5*3600) & (t<t_shut_off+ 1.5*3600+dt))./1e6)

load('CAESCav_Charging_L_P_7-4MPa_m_100.mat')
plot(x_f./L,P_f(:,t==t_shut_off)./1e6,'k')
legend('t_{shut off}','t_{shut off}+0.5h','t_{shut off}+1h','t_{shut off}+1.5h','Cavern (t=t_{shut off})')
ylabel('P [MPa]')
xlabel('Normalised length (x/L) [-]')
applystyle2plot
grid on

% Discharging
figure('color',[1 1 1])
colororder('sail')
load('CAESPipe_Discharging_L_P_7-4MPa_m_100.mat')
plot(x_n./L,P(:,(t>=t_shut_off) & (t<t_shut_off+ dt))./1e6)
hold on
plot(x_n./L,P(:,(t>=t_shut_off + 0.5*3600) & (t<t_shut_off+ 0.5*3600+dt))./1e6)
plot(x_n./L,P(:,(t>=t_shut_off + 1*3600) & (t<t_shut_off+ 1*3600+dt))./1e6)
plot(x_n./L,P(:,(t>=t_shut_off + 1.5*3600) & (t<t_shut_off+ 1.5*3600+dt))./1e6)

% load('CAESCav_Discharging_L_P_7-4MPa_m_100.mat')
% plot(x_f./L,P_f(:,t==t_shut_off)./1e6,'k')
legend('t_{shut off}','t_{shut off}+0.5h','t_{shut off}+1h','t_{shut off}+1.5h','Cavern (t=t_{shut off})')
ylabel('P [MPa]')
xlabel('Normalised length (x/L) [-]')
applystyle2plot
grid on

%% Different Depth of Discharge (DoD)
P_maxs = 7; % MPa
P_mins = 2:5; % MPa
ms = 100; % [kg/s]

Process_type = 'Charging_L';
etaX_stor_M = zeros(length(P_maxs),length(P_mins));
% Duration = zeros(length(Ps),length(ms));
% elapsedTimes = zeros(length(Ps),length(ms));

for ii=1:length(P_maxs)
    for jj = 1:length(P_mins)
        % filename = strcat('CAESPipe_Charging_L_P_',num2str(Ps(ii)),'MPa_m_',num2str(ms(jj)));
        % filename = strcat('CAESPipe_Charging_L_P_',num2str(Ps(ii)),'-',num2str(Pmins(ii)),'MPa_m_',num2str(ms(jj)));
        filename = strcat('CAESPipe_',Process_type,'_P_',num2str(P_maxs(ii)),'-',num2str(P_mins(jj)),'MPa_m_',num2str(ms))
        load(filename)
        % Ps(ii)
        % ms(jj)
        % t_shut_off
        etaX_stor_M(ii,jj) = etaX_stor;
        % Duration(ii,jj) = t_shut_off;
        % elapsedTimes(ii,jj) = elapsedTime;
        % figure('Color',[1 1 1]);
        % plot(abs(mean(error_hist)))
    end
end

figure('color',[1 1 1])
% colororder('sail')
plot(P_maxs - P_mins,etaX_stor_M)
xlabel('Depth of discharge (P_{max} - P_{min}) [MPa]')
ylabel('\eta_{disch}')

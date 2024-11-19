Ps = 6:10;
Pmins = Ps - 3;
ms = 50:50:200;

Process_type = 'Discharging_L';
% etaX_stor_M = zeros(length(Ps),length(ms));
% Duration = zeros(length(Ps),length(ms));
% elapsedTimes = zeros(length(Ps),length(ms));

for ii=1:length(Ps)
    for jj = 1:length(ms)

        filename = strcat('CAESPipe_',Process_type,'_P_',num2str(Ps(ii)),'-',num2str(Pmins(ii)),'MPa_m_',num2str(ms(jj)));
        load(filename)

        if strcmp(Process,'Charging_L')
            etaX_stor = ( XX(end) - XX(1) ) / sum(dXX);
        elseif strcmp(Process,'Discharging_L')
            etaX_stor = sum(dXX) / ( XX(end) - XX(1) );
        else
            error('Process not found or eta not implemented for the process.')
        end
        
        save(filename)

    end
end
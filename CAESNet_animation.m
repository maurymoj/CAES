% x = 1:40;
clear F
loops = 60;

% P_min = 5000000;
% P_max = 5500000;
% y_min = P_o./1e6 - 0.1;
% y_max = max(max(P_f))./1e6 + 0.1;
% y_min = min(min(v)) - 0.1;
% y_max = max(max(v)) + 0.1;
% y_min = min(min(rho_f));
% y_max = max(max(rho_f));
P_min = P_o./1e6 - 0.;
P_max = max(max(P))./1e6 + 0.1;
v_min = min(min(v)) - 0.1;
v_max = max(max(v)) + 0.1;

time = datetime("now");
vid_titl = strcat(simType,'_',string(year(time)),'_',string(month(time)),'_',string(day(time)),'_',string(hour(time)),'h',string(minute(time)),'min.avi');
vid = VideoWriter(vid_titl);
vid.FrameRate = 10;
open(vid)

F(loops) = struct('cdata',[],'colormap',[]);
figure('Color',[1 1 1])
% dt_f = 0.1; % time between frames
dt_f = 1;
dec = dt_f/dt;
t_dec = zeros(1,floor((length(t)-1)/dec+1));
% decimation = 1000;
for i=1:length(t_dec)
    t_dec(i) = t(1+dec*(i-1));
    % single plot
    % plot(x_n,P(:,1+dec*(i-1))./1e6);
    % xlabel('x [km]')
    % ylabel('P [MPa]')
    % ylim([P_min-0.1 P_max+0.1])
    % ylim([P_min max(P(:,1+dec*(i-1)))./1e6+0.001])
    % ylim([min(P(:,1+dec*(i-1)))/1e6-0.001 max(P(:,1+dec*(i-1)))/1e6+0.001])

    % plot(x_f,P_f(:,1+dec*(i-1))./1e6);
    % plot(x_f,v(:,1+dec*(i-1)));

    % plot(x_f,rho_f(:,1+dec*(i-1)));

    % double plot
    % yyaxis left
    % % plot(x_f,P_f(:,1+dec*(i-1))./1e6);
    % ylim([P_min P_max])
    % xlabel('x [km]')
    % ylabel('P [MPa]')
    % yyaxis right
    plot(x_f,v(:,1+dec*(i-1)));
    ylim([v_min v_max])
    ylabel('v [m/s]')

    ti = strcat('t = ',num2str(t_dec(i)),' s');
    title(ti);

    grid on
    drawnow
    F(i) = getframe(gcf);
    writeVideo(vid,F(i))
end
close(vid)
% fig = figure;
% movie(fig,F,1)
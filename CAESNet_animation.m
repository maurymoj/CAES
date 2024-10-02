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
P_lower_bound = P_min./1e6 - 0.;
P_upper_bound = max(max(P))./1e6 + 0.1;
P_min = min(min(P./1e6)) - 0.1;
P_max = max(max(P./1e6)) + 0.1;
v_min = min(min(v)) - 0.1;
v_max = max(max(v)) + 0.1;
T_min = min(min(T)) - 0.1;
T_max = max(max(T)) + 0.1;
rho_min = min(min(rho)) - 0.1;
rho_max = max(max(rho)) + 0.1;

time = datetime("now");
vid_titl = strcat(simType,'_',string(year(time)),'_',string(month(time)),'_',string(day(time)),'_',string(hour(time)),'h',string(minute(time)),'min.avi');
vid = VideoWriter(vid_titl);
vid.FrameRate = 10;
open(vid)

F(loops) = struct('cdata',[],'colormap',[]);
figure('Color',[1 1 1])
% dt_f = 0.1; % time between frames
dt_f = dt;
dec = dt_f/dt;
t_dec = zeros(1,floor((length(t)-1)/dec+1));
% decimation = 1000;

plot_type = 'Single';
% plot_type = 'Double';
Property = 'T';
for i=1:length(t_dec)
    t_dec(i) = t(1+dec*(i-1));
    
    if strcmp(plot_type,'Single')
        % single plot
        % plot(x_n,P(:,1+dec*(i-1))./1e6);
        % xlabel('x [km]')
        % ylabel('P [MPa]')
        % ylim([P_lower_bound-0.1 P_upper_bound+0.1])
        % ylim([P_lower_bound max(P(:,1+dec*(i-1)))./1e6+0.001])
        % ylim([min(P(:,1+dec*(i-1)))/1e6-0.001 max(P(:,1+dec*(i-1)))/1e6+0.001])
    
        % plot(x_f,P_f(:,1+dec*(i-1))./1e6);
        % plot(x_f,v(:,1+dec*(i-1)));
        % ylim([v_min v_max])
        % ylabel('v [m/s]')
        
        % plot(x_f,rho_f(:,1+dec*(i-1)).*v(:,1+dec*(i-1)).*A_h);
        % ylim([-200 400])
    
        % plot(x_n,T(:,1+dec*(i-1)));
        % ylim([T_min T_max])
        % ylabel('T [K]')

        if strcmp(Property,'P')
            plot(x_n,P(:,1+dec*(i-1))./1e6);
            xlabel('x [km]')
            ylabel('P [MPa]')
            % ylim([P_lower_bound-0.1 P_upper_bound+0.1])
            % ylim([P_lower_bound max(P(:,1+dec*(i-1)))./1e6+0.001])
            % ylim([min(P(:,1+dec*(i-1)))/1e6-0.001 max(P(:,1+dec*(i-1)))/1e6+0.001])
            % plot(x_f,P_f(:,1+dec*(i-1))./1e6);
            ylim([P_min P_max])
        elseif strcmp(Property,'v')
            plot(x_f,v(:,1+dec*(i-1)));
            ylim([v_min v_max])
            ylabel('v [m/s]')
        elseif strcmp(Property,'T')
            plot(x_n,T(:,1+dec*(i-1)));
            % hold on
            % plot(x_f,T_f(:,1+dec*(i-1)))
            ylim([T_min T_max])
            ylabel('T [K]')
        elseif strcmp(Property,'rho')
            plot(x_f,rho_f(:,1+dec*(i-1)));
            ylabel('Density [kg m-3]')
            ylim([rho_min rho_max])
        elseif strcmp(Property,'m_dot')
            plot(x_f,rho_f(:,1+dec*(i-1)).*v(:,1+dec*(i-1)).*A_h);
            ylabel('Mass flow rate [kg s-1]')
        end

        
        % ylim([-200 400])
        % plot(x_f,rho_f(:,1+dec*(i-1)));
    
        % double plot
        % yyaxis left
        % plot(x_f,P_f(:,1+dec*(i-1))./1e6);
        % ylim([P_lower_bound-0.1 P_upper_bound+0.1])
        % xlabel('x [km]')
        % ylabel('P [MPa]')
        % yyaxis right
        % plot(x_f,v(:,1+dec*(i-1)));
        % ylim([v_min v_max])
        % ylabel('v [m/s]')


    elseif strcmp(plot_type,'Double')

        % plot(x_f,rho_f(:,1+dec*(i-1)));
    
        % double plot
        yyaxis left
        if strcmp(Property,'Pv')
            plot(x_f,P_f(:,1+dec*(i-1))./1e6);
            % plot(x_n,P(:,1+dec*(i-1))./1e6);
            ylim([P_lower_bound-0.1 P_upper_bound+0.1])
            xlabel('x [km]')
            ylabel('P [MPa]')
            yyaxis right
            % plot(x_n,T(:,1+dec*(i-1)));
            % ylim([T_min T_max])
            % ylabel('T [K]')
            plot(x_f,v(:,1+dec*(i-1)));
            ylim([v_min v_max])
            ylabel('v [m/s]')
        elseif strcmp(Property,'TT')
            plot(x_f,T_f(:,1+dec*(i-1)));
            % plot(x_n,P(:,1+dec*(i-1))./1e6);
            ylim([T_min T_max])
            xlabel('x [km]')
            ylabel('T [K]')
            yyaxis right
            plot(x_n,T(:,1+dec*(i-1)));
            ylim([T_min T_max])
            ylabel('T [K]')
        end
    end

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
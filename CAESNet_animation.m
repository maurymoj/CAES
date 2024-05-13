% x = 1:40;
clear F
loops = 60;

P_min = 4000000;
P_max = 14000000;

sim = 'CAESP';
time = datetime("now")
vid_titl = strcat(sim,'_',string(year(time)),'_',string(month(time)),'_',string(day(time)),'_',string(hour(time)),'h.avi')
vid = VideoWriter(vid_titl);
vid.FrameRate = 10;
open(vid)

F(loops) = struct('cdata',[],'colormap',[])
figure('Color',[1 1 1])

for i=1:100
    plot(x,P(:,i*10));
    ti = strcat('t = ',num2str(i*1),' s');
    title(ti);
    ylim([P_min P_max])
    grid on
    drawnow
    F(i) = getframe(gcf);
    writeVideo(vid,F(i))
end
close(vid)
% fig = figure;
% movie(fig,F,1)
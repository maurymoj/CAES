% x = 1:40;
clear F
loops = 60;

vid = VideoWriter("Pt.avi");
vid.FrameRate = 10;
open(vid)

F(loops) = struct('cdata',[],'colormap',[])
figure('Color',[1 1 1])

for i=1:360
    plot(x,P(:,i*100));
    ti = strcat('t = ',num2str(i*10),' s');
    title(ti);
    ylim([5000000 5400000])
    grid on
    drawnow
    F(i) = getframe(gcf);
    writeVideo(vid,F(i))
end
close(vid)
% fig = figure;
% movie(fig,F,1)
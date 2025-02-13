%% Convergence speed monitoring
% [count==max_iter slow_conv good_conv fast_conv]
whos conv_speed; [sum(count==max_iter) sum(conv_speed < 0) sum(conv_speed == 0) sum(conv_speed > 0)]

%% Figures of 

figure('Color',[1 1 1]);
plot(abs(max(error_hist)))
xlabel('Iteration')
ylabel('Absolute max of the error')
grid on

figure('Color',[1 1 1]);
plot(abs(max(error_hist_v)))
xlabel('Iteration')
ylabel('Absolute max of the error')
grid on
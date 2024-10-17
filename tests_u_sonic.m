tic 

for i=1:size(P,1)
    u_sonic_Coolprop(i) = CP.PropsSI('speed_of_sound','P',P_f(i,j),'D',rho_f(i,j),'Air');
end

toc

tic
for i=1:size(P,1)
    % u_sonic_EoS(i) = sqrt(1.4*P_f(i,j)/rho_f(i,j));
    u_sonic_EoS(i) = sqrt(1.4*8.31*T_f(i,j)/0.02897);
end
toc

figure;
plot(u_sonic_EoS,u_sonic_Coolprop)
plot((u_sonic_EoS-u_sonic_Coolprop)./u_sonic_EoS)

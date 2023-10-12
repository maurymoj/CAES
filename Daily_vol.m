load('Daily_volume.mat')

d = ApplicableFor(1:end-1)-ApplicableFor(2:end);
day = ApplicableFor(d~=0);
V_d_Shell = Value(d~=0);

d1 = ApplicableFor1(1:end-1)-ApplicableFor1(2:end);
%day1 = ApplicableFor1(d1~=0);
V_d_Mobil = Value1(d1~=0);

d2 = ApplicableFor2(1:end-1)-ApplicableFor2(2:end);
%day2 = ApplicableFor2(d2~=0);
V_d_NSMDP = Value2(d2~=0);

% Millions of cubic meters

plot(day,V_d_Mobil,'DisplayName','V_d_Mobil');
hold on;
plot(day,V_d_NSMDP,'DisplayName','V_d_NSMDP');
plot(day,V_d_Shell,'DisplayName','V_d_Shell');
V_total = V_d_Mobil+V_d_NSMDP+V_d_Shell;
plot(day,V_total)
hold off;
legend('V_{Mobil}','V_{NSMDP}','V_{Shell}','V_{total}')
grid on
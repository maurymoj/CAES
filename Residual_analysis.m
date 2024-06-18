%% Mass
dm_x = zeros(n_n,n_t);
dm_t = zeros(n_n,n_t);
res = zeros(n_n,n_t);
for j=2:n_t
    for i=1:n_n 
        dm_x(i,j) = ( rho_f(i,j)*v(i,j) - rho_f(i+1,j)*v(i+1,j) )*A_h*dt;
        dm_t(i,j) = ( rho(i,j) - rho(i,j-1) )*A_h*dx;
        res(i,j) = (dm_t(i,j) - dm_x(i,j));
    end
end

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(x_n,res(:,floor(n_t/5)))
plot(x_n,res(:,floor(2*n_t/5)))
plot(x_n,res(:,floor(3*n_t/5)))
plot(x_n,res(:,floor(4*n_t/5)))
plot(x_n,res(:,n_t))
legend(num2str([floor(n_t/5);floor(2*n_t/5);...
         floor(3*n_t/5);floor(4*n_t/5);n_t]))
title('Mass residual x Position for different time instants')

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(t,res(1,:))
plot(t,res(floor(n_n/5),:))
plot(t,res(floor(2*n_n/5),:))
plot(t,res(floor(3*n_n/5),:))
plot(t,res(floor(4*n_n/5),:))
plot(t,res(n_n,:))
legend('1',num2str([floor(n_n/5);floor(2*n_n/5);...
         floor(3*n_n/5);floor(4*n_n/5);n_n]))
title('Mass residual x Time for different positions')

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(t,res(1,:))
plot(t,res(2,:))
plot(t,res(3,:))
plot(t,res(4,:))
legend('1','2','3','4')
title('Mass residual for the four first nodes')
%% Energy
dE_x = zeros(n_n,n_t);
dE_t = zeros(n_n,n_t);
res_E = zeros(n_n,n_t);
for j=2:n_t
    for i=1:n_n 
        dE_x(i,j) = ( rho_f(i,j)*v(i,j)*cp(j)*T_f(i,j) ...
            - rho_f(i+1,j)*v(i+1,j)*cp(j)*T_f(i+1,j))*A_h*dt;
        dE_t(i,j) = ( rho(i,j)*cp(j)*T(i,j) ...
            - rho(i,j-1)*cp(j)*T(i,j-1))*A_h*dx;
        res_E(i,j) = (dE_t(i,j) - dE_x(i,j));
    end
end

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(x_n,res_E(:,floor(n_t/5)))
plot(x_n,res_E(:,floor(2*n_t/5)))
plot(x_n,res_E(:,floor(3*n_t/5)))
plot(x_n,res_E(:,floor(4*n_t/5)))
plot(x_n,res_E(:,n_t))
legend(num2str([floor(n_t/5);floor(2*n_t/5);...
         floor(3*n_t/5);floor(4*n_t/5);n_t]))
title('Energy residual x Position for different time instants')

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(t,res_E(1,:))
plot(t,res_E(floor(n_n/5),:))
plot(t,res_E(floor(2*n_n/5),:))
plot(t,res_E(floor(3*n_n/5),:))
plot(t,res_E(floor(4*n_n/5),:))
plot(t,res_E(n_n,:))
legend('1',num2str([floor(n_n/5);floor(2*n_n/5);...
         floor(3*n_n/5);floor(4*n_n/5);n_n]))
title('Energy residual x time for different positions')

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(t,res_E(1,:))
plot(t,res_E(2,:))
plot(t,res_E(3,:))
plot(t,res_E(4,:))
legend('1','2','3','4')
title('Energy residual x time for four first nodes')

figure('Color',[1 1 1])
% plot(x_n,res(:,1))
hold all
plot(t(floor(n_t/3):end),res_E(1,(floor(n_t/3):end)))
plot(t(floor(n_t/3):end),res_E(floor(n_n/5),(floor(n_t/3):end)))
plot(t(floor(n_t/3):end),res_E(floor(2*n_n/5),floor(n_t/3):end))
plot(t(floor(n_t/3):end),res_E(floor(3*n_n/5),floor(n_t/3):end))
plot(t(floor(n_t/3):end),res_E(floor(4*n_n/5),floor(n_t/3):end))
plot(t(floor(n_t/3):end),res_E(n_n,floor(n_t/3):end))
legend('1',num2str([floor(n_n/5);floor(2*n_n/5);...
         floor(3*n_n/5);floor(4*n_n/5);n_n]))
title('Energy residual x time for different positions and t > dt/3')
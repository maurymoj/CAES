i = 1;

max_iter = 10;
min_iter = 4;
P = 1:10;
V = 10:-1:1;
err = 10;
tol = 2;
while (i < max_iter && err > tol) ...
            || (i <= min_iter)

    err = V(i) - P(i);
    [i err]
    i = i + 1;
    
end
%%
% clc
clear

v = 1:10;
v=v;
v_n = 400*ones(1,length(v)-1);

% v_n(:) = (v(1:end-1) >= 0).*v((1:end-1)) ...
%             +      (v(1:end-1) <  0).*v((2:end));

v_n(1) = v(2);

v_n(2:end) = (v(2:end-1) >= 0).*v((2:end-1)) ...
            +      (v(2:end-1) <  0).*v((3:end));

{v;v_n}
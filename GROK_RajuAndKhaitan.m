

% MATLAB code for simulating CAES cavern dynamics based on Raju and Khaitan (2012)
% This script simulates the discharge test to validate the heat transfer coefficients
% and can be adapted for daily operation.

clear; clc; close all;

% Constants
V = 300000;          % Cavern volume [m^3]
R = 287;             % Specific gas constant for air [J/(kg·K)]
cp = 1005;           % Specific heat at constant pressure [J/(kg·K)]
cv = cp - R;         % Specific heat at constant volume [J/(kg·K)]
Tamb = 50 + 273;     % Cavern wall temperature [K]
a = 0.2356;          % Constant for heff [W/(m^3·K)]
b = 0.0149;          % Coefficient for heff [W/(m^3·K) / (kg/s)^0.8]

% For discharge test:
% Initial conditions
T0 = 33 + 273;       % Initial temperature [K]
P0 = 70 * 1e5;       % Initial pressure [Pa] (approx 70 bar)
rho0 = P0 / (R * T0);% Initial density [kg/m^3]

% Time span for discharge test: 0 to 16 hours, in seconds
t_start = 0;
t_end = 16 * 3600;   % [s]
tspan = [t_start, t_end];

% Initial state vector y = [T; P]
y0 = [T0; P0];

% ODE solver
[t, y] = ode45(@cavern_dynamics, tspan, y0);

% Convert time to hours
t_hr = t / 3600;

% Plot results
figure;
subplot(2,1,1);
plot(t_hr, y(:,2)/1e5, 'b-', 'LineWidth', 1.5);
xlabel('Time (hr)');
ylabel('Pressure (bar)');
title('Cavern Pressure vs Time');
grid on;

subplot(2,1,2);
plot(t_hr, y(:,1) - 273, 'r-', 'LineWidth', 1.5);
xlabel('Time (hr)');
ylabel('Temperature (°C)');
title('Cavern Temperature vs Time');
grid on;

% ODE function
function dydt = cavern_dynamics(t, y)
    % y(1) = T [K]
    % y(2) = P [Pa]
    R = 287;
    a = 0.2356;          % Constant for heff [W/(m^3·K)]
    b = 0.0149;          % Coefficient for heff [W/(m^3·K) / (kg/s)^0.8]
    V = 3e5;
    cp = 1005;
    Tamb = 50 + 273;
    cv = cp - R;         % Specific heat at constant volume [J/(kg·K)]

    % Get mass flows (user-defined function, adapt as needed)
    [min_flow, mout_flow, Tin] = get_mass_flows(t / 3600);  % t in hours
    
    m_net = min_flow - mout_flow;  % Net mass flow rate [kg/s]
    
    % Density
    rho = y(2) / (R * y(1));
    
    % Effective heat transfer coefficient
    abs_m = abs(m_net);
    heff = a + b * abs_m^0.8;
    
    % d rho / dt = m_net / V
    drho_dt = m_net / V;
    
    % dT/dt from energy balance
    term1 = - (min_flow / V) * cp * (y(1) - Tin);
    term2 = R * y(1) * drho_dt;
    term3 = - heff * (y(1) - Tamb);
    dT_dt = (term1 + term2 + term3) / (rho * cv);
    
    % dP/dt = R * T * drho_dt + (P / T) * dT_dt
    dP_dt = R * y(1) * drho_dt + (y(2) / y(1)) * dT_dt;
    
    dydt = [dT_dt; dP_dt];
end

% User-defined mass flow function (for discharge test)
% min_flow: incoming mass flow [kg/s]
% mout_flow: outgoing mass flow [kg/s]
% Tin: incoming air temperature [K]
function [min_flow, mout_flow, Tin] = get_mass_flows(t_hr)
    % For discharge test (venting to atmosphere)
    min_flow = 0;  % No incoming during discharge
    Tin = 50 + 273;  % Default, not used during discharge
    
    if t_hr <= 5
        mout_flow = 417;
    elseif t_hr <= 15
        % Approximate linear decrease from 417 to 0 between 5-15 hr
        mout_flow = 417 * (1 - (t_hr - 5) / 10);
    else
        mout_flow = 0;
    end
end

% Note: To simulate the daily schedule, modify get_mass_flows(t_hr) based on Fig.7
% For example, define piecewise mass flows for compressor (positive m_net) and turbine (negative).
% During charging, min_flow = positive value, mout_flow=0, Tin=50+273
% During discharging, min_flow=0, mout_flow=positive, Tin irrelevant

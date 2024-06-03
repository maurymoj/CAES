function tank_filling_simulation()
    % Constants
    tank_volume = 44000;  % cubic meters
    daily_flow_rate = 17e6;  % standard cubic meters of air per day
    pipeline_length = 1000;  % meters
    pipeline_diameter = 0.9;  % meters
    initial_pressure = 0;  % ambient pressure in Pa
    final_pressure = 7e6;  % pressure in Pa

    % Simulation parameters
    time_span = [0, 1];  % simulation time span in days
    initial_conditions = [initial_pressure, 0];  % initial conditions [pressure, volume]

    % Solving the differential equation
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [time, solution] = ode45(@tank_filling_ode, time_span, initial_conditions, options);

    % Extracting results
    pressure = solution(:, 1);
    volume = solution(:, 2);

    % Plotting results
    figure;
    subplot(2, 1, 1);
    plot(time, pressure / 1e6);  % Convert pressure to MPa for plotting
    title('Pressure in the Tank Over Time');
    xlabel('Time (days)');
    ylabel('Pressure (MPa)');
    grid on;

    subplot(2, 1, 2);
    plot(time, volume);
    title('Volume in the Tank Over Time');
    xlabel('Time (days)');
    ylabel('Volume (cubic meters)');
    grid on;

    % Display final volume and pressure
    fprintf('Final volume in the tank: %.2f cubic meters\n', volume(end));
    fprintf('Final pressure in the tank: %.2f MPa\n', pressure(end) / 1e6);
end

function dydt = tank_filling_ode(t, y)
    % ODE function describing the tank filling process
    % y(1): pressure, y(2): volume

    % Constants
    daily_flow_rate = 17e6;  % standard cubic meters of air per day
    pipeline_length = 1000;  % meters
    pipeline_diameter = 0.9;  % meters

    % Calculating flow rate into the tank
    flow_rate = pipeline_flow_rate(y(1), pipeline_length, pipeline_diameter);

    % ODEs
    dydt = zeros(2, 1);
    dydt(1) = flow_rate / y(2);  % dp/dt = Q/V
    dydt(2) = daily_flow_rate;  % dV/dt = Q (constant daily flow rate)
end

function q = pipeline_flow_rate(p, L, d)
    % Function to calculate the flow rate in the pipeline using the Hagen–Poiseuille equation
    % p: pressure, L: length, d: diameter

    % Constants
    viscosity = 1.846e-5;  % dynamic viscosity of air at room temperature in Pa s

    % Hagen–Poiseuille equation
    q = (pi * (d^4) / 128 / viscosity / L) * (p - 0);  % Assuming initial pressure is 0
end
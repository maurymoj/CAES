    % CAES Simulation: Models charging, discharging, or cycling of air in a cylindrical cavern.
    % Assumptions: Ideal gas, isothermal or adiabatic conditions, constant mass flow rates.
    % Inputs: Process ('Charging', 'Cycle', 'Discharging'), heat_transfer_model ('Isothermal', 'Adiabatic')
    % Outputs: Pressure (P), Temperature (T), Exergy (X), and plots.
    clear
    clc
    
    CP = py.importlib.import_module('CoolProp.CoolProp');
    
    %--------------------------- Setup ---------------------------------
    
    % CAES process
    % Options:
    % 'Charging';
    % 'Cycle';
    % 'Idle'
    % 'Discharging';
    Process = 'Charging';
    
    % heat_transfer_model 
    % 'Adiabatic'
    % 'Isothermal'
    % Not fully implemented yet:
    % 'Steady_state'
    
    heat_transfer_model = ['Isothermal'];
    
    %----------------------- PROBLEM PARAMETERS -------------------------%
    if strcmp(Process,'Charging') || strcmp(Process,'Discharging')
        Dt = 8*3600;
    elseif strcmp(Process,'Cycle')
        Dt_charg = 8*3600;
        Dt = 16*3600;
    else
        error('Process not identified.')
    end
    
        % Cavern dimensions (assumed cylindrical)
    H = 50.625; % Height - To achieve volume similar to 0.9m, 100 km pipeline
    D = 40; % Diameter
    
    % Ambient conditions
    P_a = 101325;
    T_a = 273.15 + 25;
    rho_a = CP.PropsSI('D','P',P_a,'T',T_a,'Air');
    
    % A - Left side - Inlet 
    % P_max = 7e6;
    % DoD = 3e6; % Depth of discharge (in terms of pressure in Pa)
    % P_min = P_max - DoD;
    
    % P_max = 6.95941e6;
    % P_min = 4.056e6;
    P_max = 7e6;
    P_min = 4e6;
    
    if strcmp(Process,'Charging') || strcmp(Process,'Cycle') 
        P_in = P_max;
        T_in = 273.15 + 60;
    elseif strcmp(Process,'Discharging') 
        
    else
        error('Process not identified')
    end
    
    % m_in = 108; % Huntorf
    % m_in = rho_a*Q_a;
    m_in = 100;
    m_out = 100; % absolute value, sign is added later
    
    % Ambient conditions
    P_amb = 101325;
    T_amb = 273.15 + 25;
    
    % Dead state
    P_o = P_amb;
    T_o = T_amb;
    
    R = 287;
    g = 9.81;
    theta = 0;
    
    %--------------------- SIMULATION PARAMETERS ------------------------%
    dt = 10; % [s]
    
    
    %---------------------- ARRAYS INITIALIZATION ----------------------%
    t = 0:dt:Dt;
    
    n_t = length(t);  % n of time steps
    
    P = zeros(n_t,1);
    T = zeros(n_t,1);
    rho = zeros(n_t,1);
    cp = zeros(n_t,1);
    v = zeros(n_t,1);
    u = zeros(n_t,1);
    s = zeros(n_t,1);
    h_flow = zeros(n_t,1);
    
    % Sanity check
    m = zeros(n_t,1); % Total mass in pipeline
    E = zeros(n_t,1); % Total energy in pipeline
    
    if strcmp(Process,'Charging') || strcmp(Process,'Cycle')
        P_0 = P_min;
        T_0 = T_amb;
        % T_0 = 278.15;
        m_dot = m_in;
        stage = 'Charging';
    elseif strcmp(Process,'Discharging')
        P_0 = P_max;
        T_0 = T_amb;
        m_dot = -m_out;
        stage = 'Discharging';
    else
        error('A valid process was not selected. Select either Charging, Discharging or Cycle.')
    end
    
    %--------------------- Properties at t=0 --------------------------%
    
    
    
    if strcmp(Process,'Charging') || strcmp(Process,'Cycle')
        rho_in = CP.PropsSI('D','P',P_in,'T',T_in,'Air');
        cp_flow = CP.PropsSI('C','P',P_in,'T',T_in,'Air');
        h_in = CP.PropsSI('H','P',P_in,'T',T_in,'Air');
    elseif strcmp(Process,'Discharging')
        
    end
    
    P(1) = P_0;
    T(1) = T_0;
    rho(1) = CP.PropsSI('D','P',P(1),'T',T(1),'Air');
    cp(1) = CP.PropsSI('C','P',P(1),'T',T(1),'Air');
    
    v(1) = 1/rho(1);
    u(1) = CP.PropsSI('U','D',rho(1),'T',T(1),'Air');
    s(1) = CP.PropsSI('S','D',rho(1),'T',T(1),'Air');
    
    v_o = 1/CP.PropsSI('D','P',P_o,'T',T_o,'Air');
    u_o = CP.PropsSI('U','P',P_o,'T',T_o,'Air');
    s_o = CP.PropsSI('S','P',P_o,'T',T_o,'Air');
    
    A_h = pi*D^2/4;
    Vol = A_h*H;
    
    m(1) = rho(1)*Vol;
    E(1) = m(1)*u(1);
    
    stage_hist = {stage};
    
    for i = 2:length(t)
        if strcmp(Process,'Charging') 
            if strcmp(stage,'Charging') & P(i-1) >= P_max
                m_dot = 0;
                i_charg_end = i;
                stage = 'Idle_charg';
                disp('Charging completed')
            end
        elseif strcmp(Process,'Cycle')
            if strcmp(stage,'Charging') & P(i-1) >= P_max
                m_dot = 0;
                i_charg_end = i;
                stage = 'Idle_charg';
                disp('Charging completed')
            elseif strcmp(stage,'Idle_charg') & t(i) >= Dt_charg
                m_dot = -m_out;
                stage = 'Discharging';
                disp('Discharging started')
            elseif strcmp(stage,'Discharging') & P(i-1) <= P_min
                m_dot = 0;
                i_disch_end = i;
                stage = 'Idle_disch';
                disp('Discharging completed')
            end
            stage_hist = [stage_hist;stage];
            % stage_hist{end+1} = stage;
        elseif strcmp(Process,'Discharging')
            if strcmp(stage,'Discharging') & P(i-1) <= P_min
                m_dot = 0;
                i_disch_end = i;
                stage = 'Idle_disch';
                disp('Discharging completed')
            end
        else
            error('Process not identified/implemented')
        end
        
        m(i) = m(i-1) + m_dot*dt;
        rho(i) = m(i)/Vol;
        
        if strcmp(stage,'Charging')
            % cp_flow = CP.PropsSI('C','P',P_in,'T',T_in,'Air');
            % E(i) = E(i-1) + m_dot*cp_flow*T_in*dt;

            h_flow(i) = h_in;
            E(i) = E(i-1) + m_dot*h_flow(i)*dt;
        elseif strcmp(stage,'Discharging')
            % cp_flow = CP.PropsSI('C','P',P(i-1),'T',T(i-1),'Air');
            % E(i) = E(i-1) + m_dot*cp_flow*T(i-1)*dt;

            h_flow(i) = CP.PropsSI('H','P',P(i-1),'T',T(i-1),'Air');
            E(i) = E(i-1) + m_dot*h_flow(i)*dt;
        elseif strcmp(stage,'Idle_charg') || strcmp(stage,'Idle_disch')
            E(i) = E(i-1);
        else
            error('Stage not identified')
        end
    
        if m(i) < 0
            error('Mass became negative at step %d', i);
        end
        
        if strcmp(heat_transfer_model,'Isothermal')
            T(i) = T_0;
        elseif strcmp(heat_transfer_model,'Adiabatic')
            u(i) = E(i) / m(i); % Specific internal energy
            T(i) = CP.PropsSI('T', 'D', rho(i), 'U', u(i), 'Air');
            % T(i) = E(i)/(m(i)*cp(i-1));
        else
            error('Heat transfer model not identified')
        end
        
        P(i) = CP.PropsSI('P','D',rho(i),'T',T(i),'Air');
        cp(i) = CP.PropsSI('C','D',rho(i),'T',T(i),'Air');
    
        v(i) = 1/rho(i);
        % u(i) = CP.PropsSI('U','P',P(i),'D',rho(i),'Air');
        u(i) = CP.PropsSI('U','D',rho(i),'T',T(i),'Air');
        s(i) = CP.PropsSI('S','D',rho(i),'T',T(i),'Air');    
        
    end
    
    if strcmp(heat_transfer_model,'Isothermal')
        X = ( m.*R*T_o.*(P_o./P-1+log(P./P_o)) )./(3600*1e6); % [MWh]
        X_coolprop = m.*( (u-u_o) + P_o*(v-v_o) -T_o*(s-s_o) )./(3600*1e6); % [MWh]
    elseif strcmp(heat_transfer_model,'Adiabatic')
        X_coolprop =  m.*( (u-u_o) + P_o*(v-v_o) -T_o*(s-s_o) )./(3600*1e6); % [MWh]
        X_ideal_gas = m.*( cp.*(T - T_o - T_o*log(T./T_o))...
            + R*T.*(P_o./P - 1) ...
            +R*T_o*log(P/P_o) )./(3600*1e6); % [MWh]
        X = X_coolprop; 
        % X = cp.*(T - To - To*log(T./To))  )
        % ( m.*R*T_o.*(P_o./P-1+log(P./P_o)) )./(3600*1e6); % [MWh]
    else
        error('Heat transfer model not identified.')
    end
    
    X_m3 = m.*X./Vol; % [MWh/m3]
    
    figure('Color',[1 1 1])
    yyaxis left
    plot(t./3600,P./1e6)
    ylabel('P [MPa]')
    yyaxis right
    plot(t./3600,X-X(1))
    ylabel('X [MWh]')
    xlabel('t [h]')
    grid on
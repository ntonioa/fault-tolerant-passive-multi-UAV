% Define a class named PevMpc that models the MPC controller of a PEV (Plug-in Electric Vehicle)
classdef PevMpc
    properties
        % Default parameters for the MPC controller
        N = 24 % Prediction horizon
        T = 1/3 % Time step duration

        x_max = 12 % Maximum state value
        x_min = 1 % Minimum state value
        x_ref = 8.1 % Reference state value

        P_ch_max = 5 % Maximum charging power
        P_ch_min = 1.3 % Minimum charging power
        P_dis_max = 5 % Maximum discharging power
        P_dis_min = 1.3 % Minimum discharging power
        eta_ch = 0.955 % Charging efficiency
        eta_dis = 0.955 % Discharging efficiency
        xi = 0.15*ones(24, 1) % Perturbation factors for the objective function

        sol = struct([]) % Solution structure
    end
    
    methods
        % Constructor for the PevMpc class
        function pevMpc = PevMpc(varargin)
            if(nargin == 12)
                % Set the parameters based on the input arguments
                pevMpc.N = varargin{1}; 
                pevMpc.T = varargin{2}; 
                pevMpc.x_max = varargin{3}; 
                pevMpc.x_min = varargin{4}; 
                pevMpc.x_ref = varargin{5}; 
                pevMpc.P_ch_max = varargin{6}; 
                pevMpc.P_ch_min = varargin{7}; 
                pevMpc.P_dis_max = varargin{8}; 
                pevMpc.P_dis_min = varargin{9}; 
                pevMpc.eta_ch = varargin{10}; 
                pevMpc.eta_dis = varargin{11}; 
                pevMpc.xi = varargin{12}; 
            end
        end

        % Perform one iteration of the MPC algorithm
        function pevMpc = pevMpcIter(pevMpc, x_1, i_ref, lambda, mu, nu)
            % Define optimization variables
            P_ch = optimvar('P_ch', pevMpc.N, 1); % Charging power
            P_dis = optimvar('P_dis', pevMpc.N, 1); % Discharging power
            P = optimvar('P', pevMpc.N, 1, 'LowerBound', -pevMpc.P_dis_max*ones(pevMpc.N, 1), 'UpperBound', pevMpc.P_ch_max*ones(pevMpc.N, 1)); % Net power
            delta_ch = optimvar('delta_ch', pevMpc.N, 1, 'Type', 'integer', 'LowerBound', zeros(pevMpc.N, 1), 'UpperBound', ones(pevMpc.N, 1)); % Binary variable for charging
            delta_dis = optimvar('delta_dis', pevMpc.N, 1, 'Type', 'integer', 'LowerBound', zeros(pevMpc.N, 1), 'UpperBound', ones(pevMpc.N, 1)); % Binary variable for discharging
            x = optimvar('x', pevMpc.N+1, 1, 'LowerBound', pevMpc.x_min*ones(pevMpc.N+1, 1), 'UpperBound', pevMpc.x_max*ones(pevMpc.N+1, 1)); % State of charge variable
            
            % Define constraints
            c1 = P == P_ch-P_dis; % Net power constraint
            c2 = P_ch <= delta_ch*pevMpc.P_ch_max; % Maximum charging power constraint
            c3 = P_dis <= delta_dis*pevMpc.P_dis_max; % Maximum discharging power constraint
            c4 = delta_ch+delta_dis <= ones(pevMpc.N, 1); % Only one charging or discharging mode can be active at a time
            c5 = delta_ch*pevMpc.P_ch_min <= P_ch; % Minimum charging power constraint
            c6 = delta_dis*pevMpc.P_dis_min <= P_dis; % Minimum discharging power constraint
            c7 = x(2:pevMpc.N+1) == x(1:pevMpc.N)+pevMpc.T*(pevMpc.eta_ch*P_ch-P_dis/pevMpc.eta_dis); % State update equation
            c8 = x(1) == x_1; % Initial state constraint
            c9 = x(i_ref:end) == pevMpc.x_ref; % Reference state constraint
            
            % Define the optimization problem
            prob = optimproblem;
            prob.Objective = (pevMpc.xi'+mu'-nu'+lambda')*P; % Objective function
            prob.Constraints.c1 = c1;
            prob.Constraints.c2 = c2;
            prob.Constraints.c3 = c3;
            prob.Constraints.c4 = c4;
            prob.Constraints.c5 = c5;
            prob.Constraints.c6 = c6;
            prob.Constraints.c7 = c7;
            prob.Constraints.c8 = c8;
            prob.Constraints.c9 = c9;
            
            % Solve the optimization problem
            opts = optimoptions(@intlinprog, 'Display', 'off');
            pevMpc.sol = solve(prob, 'Options', opts); % Store the solution
        end
    end
end
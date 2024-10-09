classdef PevMpc % models the MPC controller of a PEV
    properties
        % default parameters
        N = 24
        T = 1/3

        x_max = 12
        x_min = 1
        x_ref = 8.1

        P_ch_max = 5
        P_ch_min = 1.3
        P_dis_max = 0
        P_dis_min = 0
        eta_ch = 0.955
        eta_dis = 1
        xi = 0.15*ones(1, 24)

        s_up = -inf*ones(1, 24)
        s_down = inf*ones(1, 24)
        sol = struct([])
    end
    methods
        % constructor
        function pevMpc = PevMpc(varargin) %N, T, x_max, x_min, x_ref, P_ch_max, P_ch_min, P_dis_max, P_dis_min, eta_ch, eta_dis, xi
            if(nargin == 12)
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

        % mpc iteration
        function pevMpc = pevMpcIter(pevMpc, x_1, i_ref, lambda, mu, ni)
            % variables
            P_ch = optimvar('P_ch', 1, pevMpc.N);
            P_dis = optimvar('P_dis', 1, pevMpc.N);
            P = optimvar('P', 1, pevMpc.N, 'LowerBound', -pevMpc.P_dis_max, 'UpperBound', pevMpc.P_ch_max);
            delta_ch = optimvar('delta_ch', 1, pevMpc.N, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
            delta_dis = optimvar('delta_dis', 1, pevMpc.N, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
            x = optimvar('x', 1, pevMpc.N+1, 'LowerBound', pevMpc.x_min, 'UpperBound', pevMpc.x_max);
            
            % constraints
            c1 = P == P_ch-P_dis;
            c2 = P_ch <= delta_ch*pevMpc.P_ch_max;
            c3 = P_dis <= delta_dis*pevMpc.P_dis_max;
            c4 = delta_ch+delta_dis <= 1;
            c5 = delta_ch*pevMpc.P_ch_min <= P_ch;
            c6 = delta_dis*pevMpc.P_dis_min <= P_dis;
            c7 = x(2:pevMpc.N+1) == x(1:pevMpc.N)+pevMpc.T*(pevMpc.eta_ch*P_ch-P_dis/pevMpc.eta_dis);
            c8 = x(1) == x_1;
            c9 = x(i_ref:end) == pevMpc.x_ref;
            
            % problem
            prob = optimproblem;
            prob.Objective = P*(pevMpc.xi'+mu-ni+lambda);
            prob.Constraints.c1 = c1;
            prob.Constraints.c2 = c2;
            prob.Constraints.c3 = c3;
            prob.Constraints.c4 = c4;
            prob.Constraints.c5 = c5;
            prob.Constraints.c6 = c6;
            prob.Constraints.c7 = c7;
            prob.Constraints.c8 = c8;
            prob.Constraints.c9 = c9;
            
            % solution
            opts = optimoptions(@intlinprog, 'Display', 'off');
            pevMpc.sol = solve(prob, 'Options', opts);
            pevMpc.s_up = max(pevMpc.s_up, pevMpc.sol.P);
            pevMpc.s_down = min(pevMpc.s_down, pevMpc.sol.P);
        end
    end
end
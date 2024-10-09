function sol = CenMpc(N, T, m, x_min, x_max, x_1, x_ref, i_ref, P_ch_min, P_ch_max, P_dis_min, P_dis_max, P_max, P_ref, eta_ch, eta_dis, xi)
    %% variables
    x = optimvar ('x', m, N+1, 'LowerBound', repmat(x_min, 1, N+1), 'UpperBound', repmat(x_max, 1, N+1));
    P_ch = optimvar('P_ch', m, N);
    P_dis = optimvar('P_dis', m, N);
    P = optimvar('P', m, N, 'LowerBound', repmat(-P_dis_max, 1, N), 'UpperBound', repmat(P_ch_max, 1, N));
    t = optimvar('t', 1, N);
    delta_ch = optimvar('delta_ch', m, N, 'Type', 'integer', 'LowerBound', zeros(m, N), 'UpperBound', ones(m, N));
    delta_dis = optimvar('delta_dis', m, N, 'Type', 'integer', 'LowerBound', zeros(m, N), 'UpperBound', ones(m, N));

    %% constraints
    c1 = P == P_ch-P_dis;
    c2 = P_ch <= delta_ch.*repmat(P_ch_max, 1, N);
    c3 = P_dis <= delta_dis.*repmat(P_dis_max, 1, N);
    c4 = delta_dis+delta_ch <= ones(m, N);
    c5 = delta_ch.*repmat(P_ch_min, 1, N) <= P_ch;
    c6 = delta_dis.*repmat(P_dis_min, 1, N) <= P_dis;
    c7 = x(:, 2:N+1) == x(:, 1:N)+T*(repmat(eta_ch, 1, N).*P_ch(:, 1:N)-P_dis(:, 1:N)./repmat(eta_dis, 1, N));
    c8 = x(:, 1) == x_1;
    c9 = optimconstr(m, N);
    for i = 1:m
            c9(i, 1:N-i_ref(i)+2) = x(i, i_ref(i):end) == repmat(x_ref(i), 1, N-i_ref(i)+2);
    end
    c10 = -t <= sum(P, 1)-P_ref;
    c11 = sum(P, 1)-P_ref <= t;
    c12 = sum(P, 1) <= P_max;

    %% problem
    prob = optimproblem;
    prob.Objective = sum(t+sum(P.*xi, 1), 2);
    
    prob.Constraints.c1 = c1;
    prob.Constraints.c2 = c2;
    prob.Constraints.c3 = c3;
    prob.Constraints.c4 = c4;
    prob.Constraints.c5 = c5;
    prob.Constraints.c6 = c6;
    prob.Constraints.c7 = c7;
    prob.Constraints.c8 = c8;
    prob.Constraints.c9 = c9;
    prob.Constraints.c10 = c10;
    prob.Constraints.c11 = c11;
    prob.Constraints.c12 = c12;

    %% solution
    opts = optimoptions(@intlinprog, 'Display', 'off');
    sol = solve(prob);%, 'Options', opts);
end
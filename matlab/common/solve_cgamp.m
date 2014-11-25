function [mses, a, c, r, sig] = solve_cgamp( y, F, x, ...
    channel_type, channel_prmts, prior_dist, prior_prmts, ...
    init, t_max, diff )
    [m, n] = size(F);
    sqrF = abs(F) .^ 2;

    % Select prior
    prior = @prior_gb;

    % Select channel
    if strcmp(channel_type, 'pr')
        channel = @channel_pr;
        damp = channel_prmts(2);
    else
        channel = @channel_gaussian;
    end

    % Initialize variables
    a = init(1:n);
    c = init((n + 1):end);
    g = zeros(m, 1);

    % Main loop
    mses = [];
    for t = 1:t_max
        a_old = a;
        if t > 1
            g_old = g;
            dg_old = dg;
        end

        % Update {v, w}, {g, dg}
        v = sqrF * c; 
        w = F * a - v .* g;
        [g, dg, channel_prmts] = channel(y, w, v, channel_prmts);
        if t > 1
            g = damp * g_old + (1 - damp) * g;
            dg = damp * dg_old + (1 - damp) * dg;
        end

        % Update {sig, r}, {a, c}
        sig = -1 ./ ( sqrF' * dg );
        r = a + sig .* ( F' * g );

        [a, c] = prior(r, sig, prior_prmts);
        if norm(a - a_old, 1) / n < diff
            break
        end

        % Accum./print status
        mses = [mses; sum((y - abs(F * a)) .^ 2)];
    end
end

% CHANNELS
function [g, dg, prmts_new] = channel_gaussian( y, w, v, prmts )
    delta = prmts(1);
    
    g = (y - w) ./ (delta + v);
    dg = -1 ./ (delta + v);
    
    prmts_new = [delta * sum(abs(g) .^ 2) / sum(-dg)];
end

function [g, dg, prmts_new] = channel_pr( y, w, v, prmts )
    v_eff = v + prmts(1);

    phi = 2 * abs(y) .* abs(w) ./ v_eff;
    R0 = besseli(1, phi, 1) ./ besseli(0, phi, 1);

    g = (abs(y) .* exp(i * phase(w)) .* R0 - w) ./ v_eff;
    dg = y.^2 .* (1 - R0.^2) ./ v_eff.^2 - 1 ./ v_eff;

    %prmts_new = [delta * sum(abs(g) .^ 2) / sum(-dg)];
    prmts_new = prmts;
end

% PRIORS
function [a, c] = prior_gb( r, sig, prmts )
    rho = prmts(1);
    pr_mean = prmts(2);
    pr_var = prmts(3);

    isv = 1 ./ (pr_var + sig);
    rsc = .5 .* abs(pr_mean - r) .^ 2 .* isv;
    eff = (pr_mean .* sig + r .* pr_var) .* isv;
    vrp = pr_var .* sig .* isv;

    gamma = ((1. - rho) / rho) .* (pr_var ./ vrp) .* ...
        exp(-.5 * abs(r) .^ 2 ./ sig + rsc);

    a = eff ./ (1 + gamma);
    c = vrp ./ (1 + gamma) + .5 * gamma .* abs(a) .^ 2;
end

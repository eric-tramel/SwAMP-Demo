function [mses, a, c, r, sig] = solve_gamp( y, F, x, ...
                             channel_type, channel_prmts, prior_dist, prior_prmts, ...
                             t_max, diff )
    [m, n] = size(F);
    sqrF = F .* F;

    % Select prior
    if strcmp(prior_dist, '01')
        prior = @prior_binary;
    else
        prior = @prior_gb;
    end

    % Select channel
    if strcmp(channel_type, 'pm1')
        channel = @channel_pm1;
    else
        channel = @channel_gaussian;
    end

    % Initialize variables
    a = zeros(n, 1);
    c = ones(n, 1);
    g = zeros(m, 1);

    % Main loop
    mses = [];
    for t = 1:t_max
        a_old = a;

        % Update {v, w}, {g, dg}
        v = sqrF * c;
        w = F * a - v .* g;
        [g, dg, channel_prmts] = channel(y, w, v, channel_prmts);        

        % Update {sig, r}, {a, c}
        sig = -1 ./ ( sqrF' * dg );
        r = a + sig .* ( F' * g );

        [a, c] = prior(r, sig, prior_prmts);
        if norm(a - a_old, 1) / n < diff
            break
        end

        % Accum./print status
        mses = [mses; sum((a / norm(a) - x / norm(x)) .^ 2)];
    end
end

% CHANNELS
function [g, dg, prmts_new] = channel_gaussian( y, w, v, prmts )
    delta = prmts(1);
    
    g = (y - w) ./ (delta + v);
    dg = -1 ./ (delta + v);
    
    prmts_new = [delta * sum(g .^ 2) / sum(-dg)];
end

function [g, dg, prmts_new] = channel_pm1( y, w, v, prmts )
    v_eff = v + prmts(1);
    arg = y .* w ./ sqrt(2 .* v_eff);
    g = y ./ ( sqrt(.5 * pi .* v_eff) .* erfcx(-arg) );
    dg = -g .* (w ./ v_eff + g);

    %prmts_new = [delta * sum(g. ^ 2) / sum(-dg)];
    prmts_new = prmts;
end

% PRIORS
function [a, c] = prior_gb( r, sig, prmts )
    rho = prmts(1);
    pr_mean = prmts(2);
    pr_var = prmts(3);

    isv = 1 ./ (pr_var + sig);
    rsc = .5 .* (pr_mean - r) .* (pr_mean - r) .* isv;
    eff = (pr_mean .* sig + r .* pr_var) .* isv;
    vrp = pr_var .* sig .* isv;

    gamma = ((1. - rho) / rho) .* sqrt(pr_var ./ vrp) .* ...
        exp(-.5 * r .* r ./ sig + rsc);

    a = eff ./ (1 + gamma);
    c = bsxfun( @max, gamma .* a .^ 2 + vrp ./ (1 + gamma), 1e-19 );
end

function [a, c] = prior_binary( r, sig, prmts )
    rho = prmts(1);

    z = rho + (1 - rho) .* exp(.5 * (1 - 2 * r) ./ sig);
    a = rho ./ z;
    c = a .* (1 - a);
end

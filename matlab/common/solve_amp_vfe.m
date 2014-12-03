function [mses, a, c, r, sig,VFE] = solve_amp_vfe( y, F, x, ...
                             delta, learn_delta, prior_dist, prmts, ...
                             t_max, diff)
    [m, n] = size(F);
    sqrF = F .* F;

    % Select prior
    if prior_dist == '01'
        prior = @prior_binary;
    else
        prior = @prior_gb;
    end

    % Initialize variables
    a = zeros(n, 1);
    c = ones(n, 1);
    D = zeros(m, 1);
    V = ones(m, 1);

    % Main loop
    mses = [];
    for t = 1:t_max
        D = (y - F * a) + (sqrF * c) .* (D ./ V);
        V = delta + sqrF * c;
        sig = 1 ./ ( sqrF' * (1 ./ V) );
        r = a + sig .* (F' * (D ./ V));        

        a_old = a;
        [a, c, logz] = prior(r, sig, prmts);
        
        if (t > 1)
            L1 = .5 * sum( (y - F * a) .^2 ) / delta;
            L2 = .5 * sum( log(delta + sqrF * c) );
            L3 = sum(logz);
            L4 = sum( .5 * (c + (a - r).^2) ./ sig ); 
            VFE(t) = L1 + L2 - (L3 + L4);
        end
        
        mses = [mses mean((a - x).^2)];

        %if learn_delta
            %delta = delta * sum((D ./ V).^2) / sum(1 ./ V);
        %end

        if abs(a - a_old) < diff
            break
        end
    end
end

% PRIORS
function [a, c, logz] = prior_gb( r, sig, prmts )
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
    logz = log(rho * sqrt(vrp ./ pr_var)) - rsc + log1p(gamma);
end

function [a, c] = prior_binary( r, sig, prmts )
    rho = prmts(1);

    z = rho + (1 - rho) .* exp(.5 * (1 - 2 * r) ./ sig);
    a = rho ./ z;
    c = a .* (1 - a);
end

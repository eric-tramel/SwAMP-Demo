function [mses, a, c, r, sig,VFE] = solve_gamp_vfe( y, F, x, ...
                             channel_type, channel_prmts, prior_dist, prior_prmts, ...
                             t_max, diff )
    [m, n] = size(F);
    sqrF = F .* F;

    incup = 1.2;
    incdn = 0.8;
    damp = 0.0;
    damp_max = 0.9999;
    damp_min = 0.0001;
    
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
    % g = zeros(m, 1);
    g = y;

    % Main loop
    mses = [];
    for t = 1:t_max
        a_old = a;

        % Update {v, w}, {g, dg}
        v = sqrF * c;
        w = F * a - v .* g;
        [g, dg] = channel(y, w, v, channel_prmts);   

%         fprintf('Attempting to calc VFE\n');
        %%% Is it possible to check the VFE at this point?
        % Step 1. Calculate a fixed-point \omega term
        tol = 1e-6;
        maxit = 20;        
        omega_damp = 0.5;
        reg = 1e-10;
        w_vfe = w;
        v_vfe = v;
        a_proj = F*a;
%         fprintf('  Converging \\Omega\n');
        VFE(1) = Inf;
        if t > 1
            for i=1:maxit
                % Get the g and g' since these are already coded
                [g_,dg_] = channel(y,w_vfe,v_vfe,channel_prmts);
                % Convert to Moments
                xi  = v_vfe.*g_ + w_vfe;
                phi = (v_vfe).^2.*dg_ + v_vfe;
                % Use Moment-based Newton's Method to find solution to w^*                
                step  = a_proj - xi;
                scale = -phi ./ v_vfe;
                w_last = w_vfe;
                w_vfe = omega_damp.*w_vfe - (1-omega_damp) .* (step)./(scale + reg);                
                % Check for convergence
                conv_check = norm(w_vfe-w_last)/m;
                conv_check2 = mean(abs(1 - w_vfe./w_last));
%                 plot(w_vfe); axis tight;
%                 drawnow;                            
%                 fprintf('    [%d] ||Err|| : %0.2e | Conv : %0.2e | Conv2 : %0.2e \n',i,conv_check,norm(w_last-w_vfe)/m,conv_check2);
                if conv_check2 < tol
                    break;
                end
            end
            % Now estimate the VFE using these values
            [~,~,logZu] = channel(y,w_vfe,v_vfe,channel_prmts);
            [~,~,logZi] = prior(r,sig,prior_prmts);
            
            mu_terms = (1./(2.*v_vfe)).*(w_vfe - a_proj).^2;
            i_terms  = (1./(2.*sig)).*(c + (a-r).^2);
            
            VFE(t) = -sum(logZu) - sum(mu_terms) - sum(logZi) - sum(i_terms);
            
            if VFE(t) > VFE(t-1)
                damp = damp.*incup;
                damp = min(damp,damp_max);
            else
                damp = damp.*incdn;
                damp = max(damp,damp_min);
            end            
        end
        %%%%
    
        fprintf('Damp : %0.2e\n',damp);
        % Update {sig, r}, {a, c}
        sig_ = -1 ./ ( sqrF' * dg );
        if t>1
            sig = damp.*sig + (1-damp).*sig_;
        else
            sig = sig_;
        end
        
        r_ = a + sig .* ( F' * g );
        if t>1
            r = damp.*r + (1-damp).*r_;
        else
            r   = r_;
        end
        
        [a, c] = prior(r, sig, prior_prmts);
        if norm(a - a_old, 1) / n < diff
            break
        end

        % Accum./print status
        mses = [mses; sum((a / norm(a) - x / norm(x)) .^ 2)];
    end
    VFE(1) = nan;
end

% CHANNELS
function [g, dg, logz] = channel_gaussian( y, w, v, prmts )
    delta = prmts(1);
    
    g = (y - w) ./ (delta + v);
    dg = -1 ./ (delta + v);    
    
    %%% The LogZ calculation from my notes
    logz = -0.5*log(2*pi*delta + 2*pi*v) - (y-w).^2./(2*(delta+v));
end

function [g, dg] = channel_pm1( y, w, v, prmts )
    v_eff = v + prmts(1);
    arg = y .* w ./ sqrt(2 .* v_eff);
    g = y ./ ( sqrt(.5 * pi .* v_eff) .* erfcx(-arg) );
    dg = -g .* (w ./ v_eff + g);
end

% PRIORS
function [a, c,logz] = prior_gb( r, sig, prmts )
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
        
    %%% The LogZ calculation from my notes
    % logz = log(1-rho) + log(rho)  - (r.^2./sig) - (pr_mean-r).^2./(2*(pr_var+sig)) - 0.5*log(2*pi*pr_var + 2*pi*sig);
    %%% The LogZ calculation from FreeOpt/gb_partition.m
    s  = sig;
    m  = pr_mean;
    v = pr_var;
    r2 = r.*r;
    vps = v + s;
    mmr = m - r;
    mmr2 = mmr.*mmr;
    fac = rho .* sqrt(s ./ vps);
    rsc = -0.5 .* mmr2 ./ vps;    

    %% Log Partition function
    zOut = log(fac) + rsc;
    zIn = ((1-rho)./fac) .* exp(-0.5.*r2./s - rsc);
    logz = zOut + log1p(zIn);
    
end

function [a, c] = prior_binary( r, sig, prmts )
    rho = prmts(1);

    z = rho + (1 - rho) .* exp(.5 * (1 - 2 * r) ./ sig);
    a = rho ./ z;
    c = a .* (1 - a);
end

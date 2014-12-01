%% Parameters
n = 512;
rho = 1./16;
alpha = 6;
delta = 0;
gamma = 50;

fprintf(' - Parameters are: N = %d, \\rho = %.2f, \\alpha = %.2f, \\gamma = %d.\n', ...
    n, rho, alpha, gamma)

k = ceil(rho * n);
m = ceil(alpha * n);

%% Generate problem
x = zeros(n, 1);
supp = randperm(n, k);
x(supp) = randn(k, 1);
F = sparse(gamma / n + randn(m, n) / sqrt(n));
y = sign(F * x + sqrt(delta)*randn(m,1));

%% Setup algorithm
% Obs.: the 'signal' option is only being passed so that the MSE may be
% evaluated at each iteration; commenting it out won't change the final
% estimate!
outfile = tempname;

opts.channelType = 'probit';
opts.channelPrmts = [delta];
opts.learnDelta = 0;
opts.priorDistr = 'gb';
opts.priorPrmts = [rho, 0.0, 1.0];
opts.learnPrior = 0;
opts.initState = [zeros(n, 1); ones(n, 1)];
opts.maxIter = 50;
opts.prec = 1e-5;
opts.damp = 0.0;
opts.display = 1;
opts.signal = x;
opts.output = outfile;

% Extra Feature options
opts.mean_removal = 1;
opts.adaptive_damp = 0;
opts.calc_vfe = 0;
opts.no_violations = 0;
opts.site_rejection = 0;


%% Run Mean Removal
fprintf(' - Running G-SwAMP... ')
tic
    opts.mean_removal = 1;
    a_sw_mr = run_swamp(y, F, opts);
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
mse_sw_mr = out(:,2);
rss_sw_mr = out(:,3);
cnv_sw_mr = out(:,4);
fprintf('Elapsed time: %.2fs, MSE: %.2f dB.\n', elapsed, 10 * log10(mse_sw_mr(end))); 

%% Run without Mean Removal
fprintf(' - Running G-SwAMP... ')
tic
    opts.mean_removal = 0;
    a_sw = run_swamp(y, F, opts);
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
mse_sw = out(:, 2);
rss_sw = out(:,3);
cnv_sw = out(:,4);
fprintf('Elapsed time: %.2fs, MSE: %.2f dB.\n', elapsed, 10 * log10(mse_sw(end))); 


%% Plot results
figure(1); clf;
    subplot(2, 1, 1);
        hold('on');
        plot(x/norm(x), 'ko'); plot(a_sw/norm(a_sw), 'rx'); 
        hold('off');
        xlim([0, n]); xlabel('i'); ylabel('x(i)');
        legend('signal', 'swAMP estimate');
        box on;

    subplot(2, 1, 2);
        semilogy(mse_sw);
        xlabel('SwAMP iter.'); ylabel('MSE');
        box on;

figure(2); clf;
    hold on;
        plot(cnv_sw_mr,'-g',   'LineWidth',1,'DisplayName','Convergence (MR)');
        plot(mse_sw_mr,'-b',   'LineWidth',1,'DisplayName','MSE (MR)');
        plot(rss_sw_mr,'-r',   'LineWidth',1,'DisplayName','RSS (residual norm) (MR)');        
        
        plot(cnv_sw,'-.g',   'LineWidth',1,'DisplayName','Convergence');
        plot(mse_sw,'-.b',   'LineWidth',1,'DisplayName','MSE');
        plot(rss_sw,'-.r',   'LineWidth',1,'DisplayName','RSS (residual norm)');        
        
    hold off;
    xlabel('Iteration');
    set(gca,'YScale','log');    
    box on;
    axis tight;
    legend('Location','SouthWest');
 

%% Parameters
n = 512;
rho = 0.44;
alpha = 0.72;
delta = 1e-8;
gamma = 50;
learn_noise = 0;

fprintf(' - Parameters are: N = %d, \\rho = %.2f, \\alpha = %.2f, \\gamma = %d.\n', ...
    n, rho, alpha, gamma)

k = ceil(rho * n);
m = ceil(alpha * n);

%% Generate problem
x = zeros(n, 1);
supp = randperm(n, k);
x(supp) = randn(k, 1);
F = sparse(gamma / n + randn(m, n) / sqrt(n));
y = F * x + sqrt(delta)*randn(m,1);

%% Setup algorithm
% Obs.: the 'signal' option is only being passed so that the MSE may be
% evaluated at each iteration; commenting it out won't change the final
% estimate!
outfile = tempname;

opts.channelType = 'gaussian';
if learn_noise
    opts.channelPrmts = [1.0];
    opts.delta = 1.0;
else 
    opts.channelPrmts = [delta];
    opts.delta = 1.0;
end
opts.learnDelta = learn_noise;
opts.learnChannel = learn_noise;
opts.priorDistr = 'gb';
opts.priorPrmts = [rho, 0.0, 1.0];
opts.learnPrior = 0;
opts.initState = [zeros(n, 1); ones(n, 1)];
opts.maxIter = 200;
opts.prec = 1e-8;
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


%% Run with SwAMP
fprintf(' - Running SwAMP-MR... ')
tic
    opts.solver = 'amp';
    a_swamp = run_swamp(y, F, opts);
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
mse_swamp = out(:,2);
rss_swamp = out(:,4);
cnv_swamp = out(:,5);
fprintf('Elapsed time: %.2fs, MSE: %.2f dB.\n', elapsed, 10 * log10(mse_swamp(end))); 

%% Run with SwGAMP
fprintf(' - Running SwGAMP-MR... ')
tic
    opts.solver = 'gamp';
    a_swgamp = run_swamp(y, F, opts);
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
mse_swgamp = out(:,2);
rss_swgamp = out(:,3);
cnv_swgamp = out(:,4);
fprintf('Elapsed time: %.2fs, MSE: %.2f dB.\n', elapsed, 10 * log10(mse_swgamp(end))); 


%% Plot results
figure(1); clf;
    hold on;
        plot(cnv_swamp,'-g',   'LineWidth',1,'DisplayName','AMP Convergence');
        plot(mse_swamp,'-b',   'LineWidth',1,'DisplayName','AMP MSE ');
        plot(rss_swamp,'-r',   'LineWidth',1,'DisplayName','AMP RSS (residual norm)');        
        
        plot(cnv_swgamp,':g',   'LineWidth',1,'DisplayName','GAMP Convergence');
        plot(mse_swgamp,':b',   'LineWidth',1,'DisplayName','GAMP MSE');
        plot(rss_swgamp,':r',   'LineWidth',1,'DisplayName','GAMP RSS (residual norm)');        
        
    hold off;
    xlabel('Iteration');
    set(gca,'YScale','log');    
    box on;
    axis tight;
    legend('Location','EastOutside');
 
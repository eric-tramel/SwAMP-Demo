%% Parameters
n = 512;
rho = 0.44;
alpha = 0.72;
delta = 1e-8;
gamma = 0;
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
y = F * x + sqrt(delta) * randn(m,1);

%% Setup algorithm
outfile = tempname;

opts.channelType = 'gaussian';
if learn_noise
    opts.channelPrmts = [1.0];
    opts.delta = 1.0;
else 
    opts.channelPrmts = [delta];
    opts.delta = delta;
end
opts.learnDelta = learn_noise;
opts.learnChannel = learn_noise;
opts.priorDistr = 'gb';
opts.priorPrmts = [rho, 0.0, 1.0];
opts.learnPrior = 0;
opts.initState = [zeros(n, 1); ones(n, 1)];
opts.maxIter = 150;
opts.prec = 0;
opts.damp = 0.0;
opts.display = 1;
opts.signal = x;
opts.output = outfile;

% Extra Feature options
opts.mean_removal = 0;
opts.adaptive_damp = 0;
opts.calc_vfe = 1;
opts.no_violations = 0;
opts.site_rejection = 0;

%% Run with SwAMP
fprintf(' - Running SwAMP... ')
tic
    opts.solver = 'amp';
    a_swamp = run_swamp(y, F, opts);
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
mse_swamp = out(:, 2);
vfe_swamp = out(:, 6);
fprintf('Elapsed time: %.2fs, MSE: %.2e dB.\n', elapsed, mse_swamp(end)); 

%% Run with SwGAMP
fprintf(' - Running SwGAMP... ')
tic
    opts.solver = 'gamp';
    a_swgamp = run_swamp(y, F, opts);
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
mse_swgamp = out(:, 2);
vfe_swgamp = out(:, 5);
fprintf('Elapsed time: %.2fs, MSE: %.2e dB.\n', elapsed, mse_swgamp(end)); 

%% Plot results
fig = figure(1); clf;
    hold on;
        plot(mse_swamp,'-b',   'LineWidth',1,'DisplayName','AMP MSE ');
        plot(vfe_swamp,'-r',   'LineWidth',1,'DisplayName','AMP VFE ');        
        
        plot(mse_swgamp,':b',   'LineWidth',1,'DisplayName','GAMP MSE');
        plot(vfe_swgamp,':r',   'LineWidth',1,'DisplayName','GAMP VFE');        
        
    hold off;
    xlabel('Iteration');
    set(gca,'YScale','log');    
    box on;
    axis tight;
    legend('Location','EastOutside');
print(fig, '-dpng', 'test.png');

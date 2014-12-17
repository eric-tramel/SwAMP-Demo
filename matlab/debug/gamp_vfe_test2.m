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

vfe_swamp_pos = vfe_swamp - min(vfe_swamp)

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
label_options = {'Interpreter','latex','FontSize',14};
title_options = {'Interpreter','latex','FontSize',16};
titlestr = sprintf('SwAMP vs. SwGAMP for AWGN channel\n $N = %d$, $\\rho = %0.2f$,$\\alpha = %0.2f$,$\\Delta = %0.2e$,$\\gamma = %0.1f$',n,rho,alpha,delta,gamma);
% MSE Comparison Plot
fig = figure(1); clf;
    hold on;                             
        plot(mse_swgamp,'-xr',   'LineWidth',1,'DisplayName','GAMP MSE');        
        plot(mse_swamp,'-b',   'LineWidth',1,'DisplayName','AMP MSE ');  
    hold off;
    xlabel('Iteration',label_options{:});
    ylabel('MSE',label_options{:});
    set(gca,'YScale','log');    
    box on;
    axis tight;
    legend('Location','NorthEast');
    title(titlestr,title_options{:});
if exist('export_fig')
    export_fig('swamp-to-swgamp-AWGN-mse.png','-png','-nocrop','-transparent','-r150');
else 
    print(fig, '-dpng', 'swamp-to-swgamp-AWGN-mse.png');
end

% VFE Comparisons Plot
fig = figure(2); clf;
    shifted_vfe_swamp   = vfe_swamp - (min(vfe_swamp)) + 1;
    shifted_vfe_swgamp  = vfe_swgamp - (min(vfe_swgamp)) + 1;
    hold on;        
        plot(shifted_vfe_swgamp,'-xr',   'LineWidth',1,'DisplayName','GAMP VFE');                
        plot(shifted_vfe_swamp, '-b',   'LineWidth',1,'DisplayName','AMP VFE '); 
    hold off;
    xlabel('Iteration',label_options{:});
    ylabel('Positive-Shifted VFE',label_options{:});
    set(gca,'YScale','log');
    box on; axis tight;
    legend('Location','NorthEast');
    title(titlestr,title_options{:});
if exist('export_fig')
    export_fig('swamp-to-swgamp-AWGN-vfe.png','-png','-nocrop','-transparent','-r150');
else 
    print(fig, '-dpng', 'swamp-to-swgamp-AWGN-vfe.png');
end


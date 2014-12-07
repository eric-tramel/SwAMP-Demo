
    %% Parameters
    gamma = 2;
    n = 2048;
    rho = 0.44;
    alpha = 0.72;
    delta = 1e-8;

    fprintf(' - Parameters are: N = %d, \\rho = %.2f, \\alpha = %.2f, \\Delta = %.2e, \\gamma = %d.\n', ...
        n, rho, alpha, delta, gamma)

    k = ceil(rho * n);
    m = ceil(alpha * n);

    %% Generate problem
    x = zeros(n, 1);
    supp = randperm(n, k);
    x(supp) = randn(k, 1);
    F = gamma / n + randn(m, n) / sqrt(n);
    w = sqrt(delta) * randn(m, 1);
    y = F * x + w;

    %% Setup algorithm
    % Obs.: the 'signal' option is only being passed so that the MSE may be
    % evaluated at each iteration; commenting it out won't change the final
    % estimate!
    outfile = tempname;

    opts.solver = 'amp';
    opts.delta = 1.0;
    opts.learnDelta = 1;
    opts.priorDistr = 'gb';
    opts.priorPrmts = [rho, 0.0, 1.0];
    opts.learnPrior = 0;
    opts.initState = [F\y; ones(n, 1)];
    opts.maxIter = 200;
    opts.prec = 1e-8;
    opts.display = 0;
    opts.signal = x;
    opts.output = outfile;

    %% Run algorithms

    fprintf(' - Running GAMP... ')
    tic
    [mse_gamp, ~,~,~,~, gamp_VFE] = solve_gamp_vfe(y, F, x, 'gaussian', [delta], 'gb', [rho, 0.0, 1.0], 500, opts.prec);
    elapsed = toc;
    fprintf('Elapsed time: %.2fs, MSE: %.2e.\n', elapsed, mse_gamp(end));
    
    fprintf(' - Running AMP... ')
    tic
    [mse_amp, ~,~,~,~, amp_VFE] = solve_amp_vfe(y, F, x, [delta], 0,'gb', [rho, 0.0, 1.0], opts.maxIter, opts.prec);
    elapsed = toc;
    fprintf('Elapsed time: %.2fs, MSE: %.2e.\n', elapsed, mse_amp(end));
    

    figure(1); clf;     
    subplot(2,1,1);
        hold on;
            plot(mse_amp,'-+r','DisplayName','AMP MSE');
            plot(mse_gamp,'-^b','DisplayName','GAMP MSE');        
        hold off;
        legend('Location','EastOutside');
        axis tight;
        ylabel('MSE'); 
        xlabel('Iteration');
        box on;
        set(gca,'YScale','log');
    subplot(2,1,2);
        minval = min([amp_VFE(:);gamp_VFE(:)]);
        shift  = minval-1;
        hold on;
            plot(amp_VFE-shift,'-r','DisplayName','AMP VFE','LineWidth',2);
            plot(gamp_VFE-shift,'-b','DisplayName','GAMP VFE','LineWidth',2);
        hold off
        legend('Location','EastOutside');
        axis tight;
        ylabel('Shifted VFE'); 
        xlabel('Iteration');
        box on;
        set(gca,'YScale','log');
    
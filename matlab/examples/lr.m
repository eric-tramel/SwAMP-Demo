function lr(eta)
    %% Parameters
    n = 2048;
    rho = 0.2;
    alpha = 0.6;
    delta = 1e-8;

    fprintf(' - Parameters are: N = %d, \\rho = %.2f, \\alpha = %.2f, \\Delta = %.2e, \\eta = %.2f.\n', ...
        n, rho, alpha, delta, eta)

    k = ceil(rho * n);
    m = ceil(alpha * n);
    r = ceil(eta * n);

    %% Generate problem
    x = zeros(n, 1);
    supp = randperm(n, k);
    x(supp) = randn(k, 1);
    P = randn(m, r); Q = randn(r, n);
    F = P * Q / n;
    w = sqrt(delta) * randn(m, 1);
    y = F * x + w;

    %% Setup algorithm
    % Obs.: the 'signal' option is only being passed so that the MSE may be
    % evaluated at each iteration; commenting it out won't change the final
    % estimate!
    outfile = tempname;

    opts.solver = 'amp_alt';
    opts.delta = 1.0;
    opts.learnDelta = 1;
    opts.priorDistr = 'gb';
    opts.priorPrmts = [rho, 0.0, 1.0];
    opts.learnPrior = 0;
    opts.initState = [zeros(n, 1); ones(n, 1)];
    opts.maxIter = 250;
    opts.prec = 1e-13;
    opts.display = 0;
    opts.signal = x;
    opts.output = outfile;

    %% Run algorithms
    fprintf(' - Running SwAMP... ')
    tic
    a_sw = swamp(y, F, opts);
    elapsed = toc;

    out = dlmread(outfile, ';', 1, 0);
    mse_sw = out(:, 2);
    fprintf('Elapsed time: %.2fs, MSE: %.2e.\n', elapsed, mse_sw(end)); 


    fprintf(' - Running AMP... ')
    tic
    [mse_amp, a_amp] = solve_amp(y, F, x, 1e-8, 0, 'gb', [rho, 0.0, 1.0], 250, 1e-13);
    elapsed = toc;
    fprintf('Elapsed time: %.2fs, MSE: %.2e.\n', elapsed, mse_amp(end));

    %% Plot results
    fig = figure(1);
    subplot(2, 2, [1 2]);
    hold('on');
    plot(x, 'ko'); plot(a_sw, 'rx'); plot(a_amp, 'b+');
    hold('off');
    xlim([0, n]); xlabel('i'); ylabel('x(i)');
    legend('signal', 'swAMP estimate', 'AMP estimate');

    subplot(2, 2, 3);
    semilogy(mse_sw);
    xlabel('SwAMP iter.'); ylabel('MSE');

    subplot(2, 2, 4);
    semilogy(mse_amp);
    xlabel('AMP iter.');
end

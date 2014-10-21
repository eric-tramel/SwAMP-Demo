function bit(gamma)
    %% Parameters
    n = 512;
    rho = 1./8;
    alpha = 3;

    fprintf(' - Parameters are: N = %d, \\rho = %.2f, \\alpha = %.2f, \\gamma = %d.\n', ...
        n, rho, alpha, gamma)

    k = ceil(rho * n);
    m = ceil(alpha * n);

    %% Generate problem
    x = zeros(n, 1);
    supp = randperm(n, k);
    x(supp) = randn(k, 1);
    F = sparse(gamma / n + randn(m, n) / sqrt(n));
    y = sign(F * x);

    %% Setup algorithm
    % Obs.: the 'signal' option is only being passed so that the MSE may be
    % evaluated at each iteration; commenting it out won't change the final
    % estimate!
    outfile = tempname;

    opts.channelType = 'bit';
    opts.channelPrmts = [0.0];
    opts.priorDistr = 'gb';
    opts.priorPrmts = [rho, 0.0, 1.0];
    opts.learnPrior = 0;
    opts.initState = [zeros(n, 1); ones(n, 1)];
    opts.maxIter = 50;
    opts.prec = 1e-5;
    opts.damp = 0;
    opts.display = 0;
    opts.signal = x;
    opts.output = outfile;

    %% Run algorithms
    fprintf(' - Running G-SwAMP... ')
    tic
    a_sw = swgamp(y, F, opts);
    elapsed = toc;

    out = dlmread(outfile, ';', 1, 0);
    mse_sw = out(:, 2);
    fprintf('Elapsed time: %.2fs, MSE: %.2f dB.\n', elapsed, 10 * log10(mse_sw(end))); 

    fprintf(' - Running GAMP... ')
    tic
    [mse_amp, a_amp] = solve_gamp(y, F, x, 'pm1', [], 'gb', [rho, 0.0, 1.0], 50, 1e-5);
    elapsed = toc;
    fprintf('Elapsed time: %.2fs, MSE: %.2f dB.\n', elapsed, 10 * log10(mse_amp(end)));

    %% Plot results
    fig = figure(1);
    subplot(2, 2, [1 2]);
    hold('on');
    plot(x / norm(x), 'ko'); plot(a_sw / norm(a_sw), 'rx'); plot(a_amp / norm(a_amp), 'b+');
    hold('off');
    xlim([0, n]); xlabel('i'); ylabel('x(i)');
    legend('signal', 'swAMP estimate', 'AMP estimate');

    subplot(2, 2, 3);
    plot(10 * log10(mse_sw));
    xlabel('G-SwAMP iter.'); ylabel('MSE');

    subplot(2, 2, 4);
    plot(10 * log10(mse_amp));
    xlabel('GAMP iteration');
end

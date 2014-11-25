function phr
    %% Parameters
    n = 256;
    rho = 1.0;
    alpha = 3;
    delta = 1e-8;
    ch_damp = 0.99;
    pr_var = 0.1;

    k = ceil(rho * n);
    m = ceil(alpha * n);

    %% Generate problem
    x = zeros(n, 1);
    supp = randperm(n, k);
    x(supp) = sqrt(.5 * pr_var) * (randn(k, 1) + i * randn(k, 1));
    F = sparse( sqrt(.5 / n) * (randn(m, n) + i * randn(m, n)) );
    w = sqrt(.5 * delta) * (randn(m, 1) + i * randn(m, 1));
    y = abs(F * x + w);

    %% Setup algorithm
    outfile = tempname;
    init = [ones(n, 1); rho * sqrt(pr_var) * ones(n, 1)];

    opts.channelType = 'cpr';
    opts.channelPrmts = [delta, ch_damp];
    opts.learnChannel = 0;
    opts.priorDistr = 'cgb';
    opts.priorPrmts = [rho, 0.0, pr_var];
    opts.learnPrior = 0;
    opts.initState = init;
    opts.maxIter = 500;
    opts.prec = 1e-13;
    opts.damp = 0;
    opts.display = 0;
    opts.signal = x;
    opts.output = outfile;

    %% Run algorithms
    fprintf(' - Running (complex) G-SwAMP... ')
    tic
    a_sw = cswgamp(y, F, opts);
    elapsed = toc;

    out = dlmread(outfile, ';', 1, 0);
    mse_sw = out(:, 2);
    fprintf('Elapsed time: %.2fs, RSS: %.2f dB.\n', elapsed, 10 * log10(mse_sw(end))); 

    fprintf(' - Running GAMP... ')
    tic
    [mse_amp, a_amp] = solve_cgamp(y, F, x, 'pr', [delta, ch_damp], 'gb', [rho, 0.0, pr_var], init, 500, 1e-13);
    elapsed = toc;
    fprintf('Elapsed time: %.2fs, RSS: %.2f dB.\n', elapsed, 10 * log10(mse_amp(end)));

    %% Plot results
    % Signal and reconstruction differ by a phase
    a_sw = a_sw .* exp(i * (phase(x) - phase(a_sw)));
    a_amp = a_amp .* exp(i * (phase(x) - phase(a_amp)));

    fig = figure(1);
    subplot(2, 2, [1 2]);
    hold('on');
    plot(real(x), imag(x), 'ko'); plot(real(a_sw), imag(a_sw), 'rx'); plot(real(a_amp), imag(a_amp), 'b+');
    hold('off');
    xlabel('Re(x)'); ylabel('Im(x)')
    xlim([min(real(a_sw)), max(real(a_sw))]); ylim([min(imag(a_sw)), max(imag(a_sw))]);
    legend('signal', 'swAMP estimate', 'AMP estimate');

    subplot(2, 2, 3);
    plot(10 * log10(mse_sw));
    xlabel('G-SwAMP iter.'); ylabel('RSS (dB)');

    subplot(2, 2, 4);
    plot(10 * log10(mse_amp));
    xlabel('GAMP iteration');
end

%% Parameters
gamma = 20;
n = 512;
% rho = 0.44;
rho = 0.4;
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



%% This is how Schniter et al approached this problem
% % Get global mean, column mean, row mean
% A1 = A.mult(ones(n,1));
%     mu = (ones(1,m)*A1)/(m*n);
%     obj.col = A.multTr(ones(m,1))/m - conj(mu)*ones(n,1);
%     obj.col2 = abs(obj.col).^2;
%     obj.gam = A1/n;
%     obj.gam2 = abs(obj.gam).^2;

% % compute scaling constants
%     Atildefro2 = ones(1,m)*obj.A.multSq(ones(n,1)) ...
%             -2*real( obj.gam'*A1 ...
%                 + ones(1,m)*obj.A.mult(obj.col) ...
%                 + (ones(1,m)*obj.gam)*(ones(1,n)*obj.col) ) ...
%             + n*(ones(1,m)*obj.gam2) + m*(ones(1,n)*obj.col2);
%     obj.b12 = min(1, sqrt(Atildefro2/(n*(ones(1,m)*obj.gam2)))); % protect agains gam2~=0
%     obj.b21 = sqrt(Atildefro2/(m*(n+obj.b12^2)));
%     obj.b13 = sqrt(Atildefro2/(n*m));
%     obj.b31 = sqrt(Atildefro2/(m*(ones(1,n)*obj.col2+obj.b13^2)));

%% What is the "Modified" matrix?
on = ones(n,1);
om = ones(m,1);

% Calculate the average scale of the rows of F
F_ms_row = mean(sum(F.^2,2)./n);
F_ms_col = mean(sum(F.^2,1)./m);

b21 = sqrt(F_ms_col);
b12 = 1;
b31 = 1;
b13 = sqrt(F_ms_col);
g  = (1/n)*F*on;		% Get \gamma
mu = (1/m)*om'*g;	% Get Mean
ct = (1/m)*om'*(F - mu*(om*on'));
Ft = F - g*on' - om*ct;

Fb = [Ft,         b12*g,   b13*om;
	  b21*on', -b21*b12,        0;
	  b31*ct,         0, -b31*b13];


yb = [y;0;0]; 		% Adding in some zeros


% % Do we do okay matching ?
% figure(1); clf;
% 	plot(sum(Fb.^2,2)/n);
% 	title('Final Row Magnitudes');

% figure(2); clf;
% 	plot(sum(Fb.^2,1)/m);
% 	title('Final Col Magnitudes');

%% Here is a second attempt, the Schniter way
col   = (1/m)*F'*om - mu*on;
col2  = abs(col).^2;
fro   = norm(Ft,'fro')^2;
b12_s = min(1,sqrt(fro./(n*om'*g)))
b21_s = sqrt(fro/(m*(n + b12_s^2)))
b13_s = sqrt(fro/(m*n))
b31_s = sqrt(fro/(m*(on'*col2 + b13_s^2)))

Fb_s = [Ft,             b12_s*g,     b13_s*om;
	    b21_s*on', -b21_s*b12_s,            0;
	    b31_s*ct,             0, -b31_s*b13_s];
Fb_s = sparse(Fb_s);

% % Do we do okay matching ?
% figure(1); clf;
% 	plot(sum(Fb_s.^2,2)/n);
% 	title('Final Row Magnitudes');

% figure(2); clf;
% 	plot(sum(Fb_s.^2,1)/m);
% 	title('Final Col Magnitudes');

outfile = tempname;
opts.solver = 'amp';
opts.delta = delta;
opts.learnDelta = 0;
opts.priorDistr = 'gb';
opts.priorPrmts = [rho, 0.0, 1.0];
opts.learnPrior = 0;
opts.initState = [zeros(n+2, 1); ones(n+2, 1)];
opts.maxIter = 1000;
opts.prec = delta;
opts.display = 1;
opts.signal = x;
opts.output = outfile;
opts.damp = 0.0;

% Extra Feature options
opts.mean_removal   = 1;
opts.adaptive_damp  = 0;
opts.calc_vfe       = 1;
opts.no_violations  = 0;
opts.site_rejection = 0;

%% Run algorithms
fprintf(' - Running SwAMP... ')
tic
a_with_aux = swamp(yb, Fb_s, opts);
a_sw = a_with_aux(1:(end-2));
elapsed = toc;

out = dlmread(outfile, ';', 1, 0);
iterations = size(out,1);
mse_sw = out(:, 2);
delta_sw = out(:,3);
rss_sw = out(:,4);
cnv_sw = out(:,5);

fprintf('Elapsed time: %.2fs, MSE: %.2e.\n', elapsed, mse_sw(end)); 

%% Plot results
figure(1); clf;
    subplot(2, 1, 1);
        hold('on');
        plot(x, 'ko'); plot(a_sw, 'rx'); 
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
        plot(cnv_sw,'-g',   'LineWidth',1,'DisplayName','Convergence');
        plot(mse_sw,'-b',   'LineWidth',1,'DisplayName','MSE');
        plot(delta_sw,'-m^','LineWidth',1,'DisplayName','\Delta Estimate');
        plot(rss_sw,'-r',   'LineWidth',1,'DisplayName','RSS (residual norm)');        
        plot(delta*ones(size(delta_sw)),':k','LineWidth',1,'DisplayName','True \Delta');
    hold off;
    xlabel('Iteration');
    set(gca,'YScale','log');    
    box on;
    axis tight;
    legend('Location','SouthWest');


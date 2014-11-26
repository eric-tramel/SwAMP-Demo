%% Parameters
gamma = 50;
n = 512;
% rho = 0.44;
rho = 0.3;
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
mu = (1/(m*n))*om'*g;	% Get Mean
ct = (1/m)*om'*(F - mu*(om*on'));
Ft = F - g*on' - om*ct;

Fb = [Ft,         b12*g,   b13*om;
	  b21*on', -b21*b12,        0;
	  b31*ct,         0, -b31*b13];


yb = [y;0;0]; 		% Adding in some zeros


% Do we do okay matching ?
figure(1); clf;
	plot(sum(Fb.^2,2)/n);
	title('Final Row Magnitudes');

figure(2); clf;
	plot(sum(Fb.^2,1)/m);
	title('Final Col Magnitudes');

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

% Do we do okay matching ?
figure(1); clf;
	plot(sum(Fb_s.^2,2)/n);
	title('Final Row Magnitudes');

figure(2); clf;
	plot(sum(Fb_s.^2,1)/m);
	title('Final Col Magnitudes');

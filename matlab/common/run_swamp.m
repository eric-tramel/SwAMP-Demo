function [a, c, r, s, final_prior] = run_swamp(y,F,opts)

% Check to see if we are in mean removal mode. 
mean_removal = 0;
if isfield(opts,'mean_removal')
	mean_removal = opts.mean_removal;
end

if mean_removal
    [m,n] = size(F);
    
    if isfield(opts,'initState')
        inita = opts.initState(1:n/2);
        initc = opts.initState((n/2+1):end);
        opts.initState = [inita; 0; 0; initc; 1; 1];
    end
	
	on = ones(n,1);
	om = ones(m,1);

	g  = (1/n)*F*on;					% Get \gamma
	mu = (1/m)*om'*g;					% Get Mean
	ct = (1/m)*om'*(F - mu*(om*on'));
	Ft = F - g*on' - om*ct;

	col   = (1/m)*F'*om - mu*on;
	col2  = abs(col).^2;
	fro   = norm(Ft,'fro')^2;
	b12   = min(1,sqrt(fro./(n*om'*g)));
	b21   = sqrt(fro/(m*(n + b12^2)));
	b13   = sqrt(fro/(m*n));
	b31   = sqrt(fro/(m*(on'*col2 + b13^2)));

	Fb    = [Ft,        b12*g,   b13*om;
		    b21*on', -b21*b12,        0;
		    b31*ct,         0, -b31*b13];

	if issparse(F)
		Fb = sparse(Fb);
	end

	F = Fb;
	y = [y; 0; 0];
end

% Call the proper SwAMP binary according to the channel being used
executable_mode = 'swamp'; 		% 0: SwAMP
							% 1: SwGAMP
							% 2: C-SwGAMP
switch opts.channelType
	case 'gaussian'
		executable_mode = 'swamp';
	case 'cpr'
		executable_mode = 'cswgamp';
	case 'cgaussian'
		executable_mode = 'cswgamp';
	otherwise
		executable_mode = 'swgamp';
end

if isfield(opts,'solver')
	switch opts.solver
		case 'amp'
			executable_mode = 'swamp';
		case 'amp_alt'		
			executable_mode = 'swamp';
		case 'gamp'
			executable_mode = 'swgamp';
		case 'cgamp'
			executable_mode = 'cswgamp';
	end
end

switch executable_mode
	case 'swamp'
		[a,c,r,s,final_prior] = swamp(y,F,opts);
	case 'swgamp'
		[a,c,r,s,final_prior] = swgamp(y,F,opts);
	case 'cswgamp'
		[a,c,r,s,final_prior] = cswgamp(y,F,opts);
end

if mean_removal
	a = a(1:(end-2));
	c = c(1:(end-2));
	r = r(1:(end-2));
	s = s(1:(end-2));
end
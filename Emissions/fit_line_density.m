function [ ffit, emgfit, f0, history, fitresults, N ] = fit_line_density( no2_x, no2_ld, varargin )
%FIT_LINE_DENSITY Fits an exponentiall modified Gaussian function to a line density
%   Expoentially modified Gaussian (EMG) functions work very well to fit
%   line densities of NO2 plumes observed from space because the two
%   components allow it to capture the plume build up or at slow winds
%   (Gaussian) and its decay downwind at faster winds (exponential decay).
%   c.f. de Foy et al., Atmos. Environ., 2014, pp. 66-77.
%
%   This function requires as input the x-coordinates in kilometers and
%   line densities (moles / km recommended, though other units *should*
%   work). A third optional input controls the level of console output
%   fmincon provides: it defaults to 'iter' meaning that fmincon will
%   output information about each step. Another option is 'final' which
%   will only output a final report.
%
%   This function will attempt to fit an EMG function to the line density
%   using fmincon. The five fitting terms (a, x_0, mu_x, sigma_x, and B) as
%   described in de Foy 2014 are returned as the first five outputs. The
%   sixth output is the EMG fit itself as a vector that will have the same
%   x-coordinates as input.
%
%   There are several additional outputs useful for debugging poor fits.
%   The seventh output (f0) is the initial guess of the five fitting
%   parameters. The eighth output (history) will contain a record of how
%   the fitting parameters were iterated during the minimization. The ninth
%   output contains additional information returned from fmincon.
%
%   Josh Laughner <joshlaugh5@gmail.com> 5 Feb 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addOptional('fmincon_output','iter',@(x) ismember(x,{'none','iter','final'}));
p.addParameter('fixed_param','',@ischar);
p.addParameter('fixed_val',[]);
p.addParameter('f0',[]);
p.parse(varargin{:});
pout=p.Results;

fmincon_output = pout.fmincon_output;
fixed_param = pout.fixed_param;
fixed_val = pout.fixed_val;
f0in = pout.f0;

% Default values
% fmincon_output = 'iter';
% fixed_param = '';
% fixed_val = [];
% 
% 
% v=1;
% while v <= numel(varargin)
%     if strcmp('fixed_param',varargin{v})
%         fixed_param = varargin{v+1};
%         v=v+2;
%     elseif strcmp('fixed_val',varargin{v})
%         fixed_val = varargin{v+1};
%         v=v+2;
%     elseif strcmp('f0',varargin{v})
%         f0in = varargin{v+1};
%         v=v+2;
%     else
%         fmincon_output = varargin{v};
%         v=v+1;
%     end
% end

E = JLLErrors;
if ~isvector(no2_x) || ~isnumeric(no2_x)
    E.badinput('no2_x must be a vector of numeric inputs')
end
if ~isvector(no2_ld) || ~isnumeric(no2_ld) || any(size(no2_ld) ~= size(no2_x))
    E.badinput('no2_ld must be a numeric vector the same size as no2_ld')
end
if ~isempty(fixed_param) && ~ismember(fixed_param,{'a','x0','mux','sx','B'})
    E.badinput('fixed_param (if given) must be one of a, x0, mux, sx, or B')
end
if ~isempty(fixed_val) && ~ischar(fixed_val) && ~isscalar(fixed_val)
    E.badinput('fixed_val must be a scalar number')
end
if xor(isempty(fixed_param),isempty(fixed_val))
    E.badinput('Both or neither of fixed_param and fixed_val must be set')
end
if ~isempty(f0in) && numel(f0in) ~= 5
    E.badinput('If an f0 is specified as input, it must have five elements')
end

% Define the fit function, which although is physically a function of x, is
% to be minimized by the variation of five parameters: a, x_0, sigma_x,
% mu_x, and B. This function is the sum of squared residuals between the
% EMG fit and the actual line density. The five variational parameters will
% be set as the five elements of the vector f, necessary for inserting them
% into the fmincon routine.

    function e = emgfxn(f,x)
%emgfxn = @(f,x) 
        e = f(1)/2 .* exp( (f(4)^2 / (2 * f(2)^2)) - (x - f(3)) ./ f(2) )...
        .* ( 1 - erf( (f(4)^2 - f(2).*(x - f(3)))./(sqrt(2) * f(4) * f(2)) ) ) + f(5);
        e(isnan(e)) = Inf;
    end
fitfxn = @(f) nansum((emgfxn(f,no2_x) - no2_ld).^2);

history.x = [];
opts = optimoptions('fmincon','Display',fmincon_output,'OutputFcn',@outfun);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best guess of initial conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(f0in)
    f0 = f0in;
else
    f0 = ones(1,5);

    % f(2) is x_0 or length scale for chemical decay. Assume a lifetime of 3 hr
    % and 5 m/s winds = 54 km for 1 lifetime
    f0(2) = 54;
    % f(1) is a or the mass of NO2 in one length scale. Tried considering this
    % the integral of line density from ~0 to one x_0 downwind, but that tends
    % to be too large by a factor of 10. 
    %      i = find(no2_x > 0, 1, 'first');
    %       j = find(no2_x > no2_x(i) + f0(2), 1, 'first');
    %       f0(1) = trapz(no2_x(i:j), no2_ld(i:j));
    % Instead, let's just try starting out as half the maximum value of the
    % line density.
    f0(1) = max(no2_ld)/2;

    % f(3) is mu_x or where the average for the Gaussian is, this should be
    % around the maximum.
    [~,m] = max(no2_ld);
    f0(3) = no2_x(m);

    % f(4) is sigma_x, the std. dev. of the gaussian. Can be defined as
    % FWHM/2.355.
    halfmax = (max(no2_ld) - no2_ld(1))/2 + no2_ld(1);
    mpre = find(no2_ld(1:m) < halfmax, 1, 'last');
    mpost = find(no2_ld > halfmax, 1, 'first');
    fwhm = abs(interp1(no2_ld(mpre:mpost), no2_x(mpre:mpost), halfmax));
    f0(4) = fwhm / 2.355;

    % f(5) is the constant term B, or background. So let's just set to be the
    % left-most line density.
    f0(5) = min(no2_ld);
end
%%%%%%%%%%%%%%%%%
% Set up bounds %
%%%%%%%%%%%%%%%%%
% Remember, f = [a, x0, mu_x, sigma_x, B]
f_lb = nan(1,5);
f_ub = nan(1,5);

% a (related to plume mass) must be > 0 to be physically meaningful
f_lb(1) = 0; f_ub(1) = Inf; %f_ub(1) = max(no2_ld)*1.5;

% x0 (length scale of chemical decay) must be > 0 to be physically
% meaningful. This term tends to cause problems in the fit when the
% algorithm thinks it sees a minimum where x0 gets very small, so I chose
% 0.5 km as the lower bound because it is ~1/10th of the resolution of the
% gridded data. This basically says that the decay must be observable at
% the resolution of the data.
f_lb(2) = 0.5; f_ub(2) = Inf;

% mu_x should be close to the position of the source, but what we can say
% for certain is that it must lie somewhere on no2_x. A stronger condition
% might be that it must lie within 1 std. dev. of the maximum of the line
% density, that may be useful if fitting continues to be foolish.
f_lb(3) = min(no2_x); f_ub(3) = max(no2_x);

% sigma_x describes the width of the Gaussian; it must be positive (and not
% just technically positive - it should have at least some width, so 0.5 km
% was chosen as ~1/10th the resolution of the gridded data) and should
% really be less than the full width of it. We will assume that the buildup
% on the left hand side represents this full width.
[~,m] = max(no2_ld);
f_lb(4) = 0.5; f_ub(4) = no2_x(m) - min(no2_x);

% B is the background value. It must be > 0 to be physically meaningful,
% and should really not exceed the maximum line density.
f_lb(5) = 0; f_ub(5) = max(no2_ld);

% Also constrain that mu_x + x0 must fall within the domain, i.e. assume
% that at least one chemical lifetime occurs within the distance studied.
A = [0 1 1 0 0];
b = max(no2_x);

if ischar(fixed_val) && strcmpi(fixed_val, 'f0')
    switch fixed_param
        case 'a'
            fixed_val = f0(1);
        case 'x0'
            fixed_val = f0(2);
        case 'mux'
            fixed_val = f0(3);
        case 'sx'
            fixed_val = f0(4);
        case 'B'
            fixed_val = f0(5);
    end
elseif ischar(fixed_val)
    E.badinput('fixed_val should be a number if not the string ''f0''')
end

switch fixed_param
    case 'a'
        inds=2:5;
        nlcon = @(f) nonlin_constr([fixed_val, f]);
        fitfxn_fix = @(f) fitfxn([fixed_val, f]);
        emgfxn_fix = @(f,x) emgfxn([fixed_val, f],x);
    case 'x0'
        inds=[1,3,4,5];
        nlcon = @(f) nonlin_constr([f(1), fixed_val, f(2:4)]);
        fitfxn_fix = @(f) fitfxn([f(1), fixed_val, f(2:4)]);
        emgfxn_fix = @(f,x) emgfxn([f(1), fixed_val, f(2:4)],x);
    case 'mux'
        inds=[1,2,4,5];
        nlcon = @(f) nonlin_constr([f(1:2), fixed_val, f(3:4)]);
        fitfxn_fix = @(f) fitfxn([f(1:2), fixed_val, f(3:4)]);
        emgfxn_fix = @(f,x) emgfxn([f(1:2), fixed_val, f(3:4)],x);
    case 'sx'
        inds=[1,2,3,5];
        nlcon = @(f) nonlin_constr([f(1:3), fixed_val, f(4)]);
        fitfxn_fix = @(f) fitfxn([f(1:3), fixed_val, f(4)]);
        emgfxn_fix = @(f,x) emgfxn([f(1:3), fixed_val, f(4)],x);
    case 'B'
        inds=1:4;
        nlcon = @(f) nonlin_constr([f, fixed_val]);
        fitfxn_fix = @(f) fitfxn([f, fixed_val]);
        emgfxn_fix = @(f,x) emgfxn([f, fixed_val],x);
    otherwise
        inds = 1:5;
        nlcon = @nonlin_constr;
        fitfxn_fix = @(f) fitfxn(f);
        emgfxn_fix = @(f,x) emgfxn(f,x);
end
f0 = f0(inds);
f_lb = f_lb(inds);
f_ub = f_ub(inds);
A = A(inds);

[fitparams, fitresults.fval, fitresults.exitFlag, fitresults.output, fitresults.lambda, fitresults.grad, fitresults.Hessian]...
    = fmincon(fitfxn_fix, f0, A, b, [], [], f_lb, f_ub, nlcon, opts);

switch fixed_param
    case 'a'
        ffinal = [fixed_val, fitparams];
    case 'x0'
        ffinal = [fitparams(1), fixed_val, fitparams(2:4)];
    case 'mux'
        ffinal = [fitparams(1:2), fixed_val, fitparams(3:4)];
    case 'sx'
        ffinal = [fitparams(1:3), fixed_val, fitparams(4)];
    case 'B'
        ffinal = [fitparams, fixed_val];
    otherwise 
        ffinal = fitparams;
end
ffit.a = ffinal(1);
ffit.x_0 = ffinal(2);
ffit.mu_x = ffinal(3);
ffit.sigma_x = ffinal(4);
ffit.B = ffinal(5);
emgfit = emgfxn_fix(fitparams, no2_x);

try
    opts = statset('derivstep',eps^(1/4).* f0);
    [N.beta, N.R, N.J, N.CovB, N.MSE, N.ErrorModelInfo] = nlinfit(no2_x,no2_ld,emgfxn_fix,f0);
    %[N.beta, N.R, N.J, N.CovB, N.MSE, N.ErrorModelInfo] = nlinfit(no2_x,no2_ld,emgfxn_fix,fitparams);
    N.emg = emgfxn_fix(N.beta,no2_x);
catch err
    if strcmp(err.identifier,'stats:nlinfit:NonFiniteFunOutput')
        N = struct;
        fprintf('nlinfit found NaNs or Inf in the function for inputs: \n\tfixed_param = %s \n\tfixed_val = %.3g \n\tf0in = %s \n\tf0 = %s\n', fixed_param, fixed_val, mat2str(f0in), mat2str(f0));
    else
        rethrow(err)
    end
end

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end

    function [c,ceq] = nonlin_constr(f)
        % Nonlinear constraint requires that the value of the exponential
        % term not exceed 20.  Since the error function is bound to [-1, 1]
        % and 1-erf is always on [0 2], this will hopefully keep the
        % function from going to NaN which happens when the exponetial term
        % goes to infinity while the 1-erf term goes to 0. This constraint
        % is purely to force fmincon to behave, the only possible physical
        % justification is that the a term should control the amount of
        % mass in the plume and not the exponetial itself.
        c = max((f(4)^2 / (2 * f(2)^2)) - (no2_x - f(3)) ./ f(2) - log(20));
        ceq = [];
    end
end


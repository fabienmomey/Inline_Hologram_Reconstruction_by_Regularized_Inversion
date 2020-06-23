function [xopt,varargout] = algoRI(x0,crit,options)
% [xopt,varargout] = algoRI(x0,crit,options)
%
%   This function performs the FISTA algorithm [1] for phase retrieval from 
%   intensity measurements from a transmittance plane zA propagated to a 
%   plane zB.
%
%   X0: initial image guess => manipulated (for the default case) as 
%       2-component image ([width,height,2]), each one corresponding 
%       respectively to the real and imaginary part of the complex 
%       deviation from the unit transmittance plane :
%
%                   T = 1 + X(:,:,1) + i X(:,:,2)
%
%   CRIT: function handle to the chosen criterion. It embeds:
%       Y: data image (square root of intensity measurements).
%       Hz: function handle to perform the propagation operator  (see 
%           getFresnelPropagation and propagationOperator functions).
%       H_z: function handle to perform the backpropagation operator (see 
%            getFresnelPropagation and propagationOperator functions).
%
%   OPTIONS: (structure of algorithm parameters)
%       * Lip:      Lipschitz constant for FISTA algorithm. If not set, 
%                   a backtracking rule will be applied [1].
%       * TYPE_OBJ: 
%           > 'dephasing': purely and weakly dephasing object.
%               => the unkown image X is real and represents the imaginary
%               part of the transmittance, which is approximately
%               considered as the targeted phase-shift image.
%               => the dimensions of X are just [width,height].
%           > 'absorbing': purely absorbing object.
%               => the unkown image X is real and represents the opposite
%               of the opacity : X = T - 1.
%               => the dimensions of X are just [width,height].
%           > 'unknown' (default): mix object.
%               => the unkown image X represents the complex 
%               deviation from the unit transmittance plane.
%               => the dimensions of X are just [width,height,2].
%       * FLAG_FIENUP:  a flag (default: false) to select de the Fienup 
%                       criterion in the case TYPE_OBJ = 'unknown'. If
%                       false, the weighted least squares criterion on 
%                       intensity [2] is selected.
%       * W:    diagonal elements of the inverse noise covariance matrix C^{-1}
%               => under hypothesis of uncorrelated noise [2]. 
%       * MU         :  hyperparameter value (default: 0.0) for 
%                       soft-tresholding (prox L1).
%                       => if TYPE_OBJ is set to 'unkown', this
%                       soft-thresholding will be applied to both the real
%                       and imaginary parts separately.
%                       This corresponds to solving the following inverse 
%                       problem with L1 regularization : 
%                             f(x) = || |1+H.x| - y ||_2^2 + mu * || x ||_1
%                       at each iteration.
%       * RCONST = [XRMIN,XRMAX]: a 2-element vector giving hard constraint
%                          parameter for the real part of X
%                           \_ TYPE_OBJ = 'dephasing'
%                               \_ default: [0,0] (X is purely imaginary)
%                           \_ TYPE_OBJ = 'absorbing'
%                               \_ default: [-1,0] (X is purely real)
%                           \_ TYPE_OBJ = 'unknown'
%                               \_ default: [-2,0] (because 0 < |T| < 1
%                                                   and -1 < cos(phi) < 1)
%       * ICONST = [XIMIN,XIMAX]: a 2-element vector giving hard constraint
%                          parameter for the imag part of X
%                           \_ TYPE_OBJ = 'dephasing'
%                               \_ default: [-1,1] (X is purely imaginary)
%                           \_ TYPE_OBJ = 'absorbing'
%                               \_ default: [0,0] (X is purely real)
%                           \_ TYPE_OBJ = 'unknown'
%                               \_ default: [-1,1] (because 0 < |T| < 1
%                                                   and -1 < sin(phi) < 1)
%       * SUPPORT    : mask of a support constraint to be enforced.
%       * MAXITER    : maximum number of iterations  (default: 100) 
%       * FLAG_COST  : compute cost at each iteration
%       * XTRUE      : ground truth. If set, output information giving the
%                      evolution of the SNR and the cost can be computed.
%       * VERBOSE    : verbose mode (default: true)
%
% References
%
% - [1]  A. Beck and M. Teboulle, “Fast gradient-based algorithms for con-
%       strained total variation image denoising and deblurring problems,”
%       IEEE Trans. Image Process. 18, 2419–2434, 2009.
% - [2] F. Momey, L. Denis, T. Olivier, C. Fournier, "From Fienup’s phase 
%                   retrieval techniques to regularized inversion for 
%                   in-line holography: tutorial," JOSA A, vol. 36, no. 12, 
%                   D62-D80, 2019. 
%
% Created: 05/27/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516,
%           Université Jean Monnet,
%           F-42000 Saint-Étienne,
%           France
%           fabien.momey@univ-st-etienne.fr
%

%% Extract size (in pixels) of the field of view
[npix_W, npix_H, ncomp] = size(x0);

%% Extract reconstruction parameters

%% Lipschitz constant L
flag_backtrack = true;
if (isfield(options, 'Lip'))
    flag_backtrack = false;
else
    options.Lip = 1.0;
    options.eta = 1.1;
end

if (flag_backtrack)
     warning('FISTA with backtracking.');
end

%% MU hyperparameter for performing a soft-thresholding 
%% (default: 0.0)
if (~isfield(options, 'mu'))
    options.mu = 0.0;
end

%% Extract hard constraints bounds
if (~isfield(options, 'rconst'))
    warning('RCONST not set: set to [-Inf,Inf] (no constraint).');
    options.rconst = [-2.0,0.0];
else
    if (any(size(options.rconst)~=[1,2]) && any(size(options.rconst)~=[2,1]))
        error('RCONST must be a 2-scalar vector.');
    end
    options.rconst=sort(options.rconst);
    options.rconst = [max(options.rconst(1),-2.0),min(options.rconst(2),0.0)];
end

if (~isfield(options, 'iconst'))
    warning('ICONST not set: set to [-Inf,Inf] (no constraint).');
    options.iconst = [-1.0,1.0];
else
    if (any(size(options.iconst)~=[1,2]) && any(size(options.iconst)~=[2,1]))
        error('ICONST must be a 2-scalar vector.');
    end
    options.iconst=sort(options.iconst);
    options.iconst = [max(options.iconst(1),-1.0),min(options.iconst(2),1.0)];
end

%% Extract TYPE_OBJ and adapt hard constraints bounds 
if (~isfield(options, 'type_obj'))
    options.type_obj = 'unknown';
end

if (strcmp(options.type_obj,'dephasing'))
    options.rconst = [0.0,0.0];
elseif (strcmp(options.type_obj,'absorbing'))
    options.iconst = [0.0,0.0];
    options.rconst = [-1.0,options.rconst(2)];
end

% %% Get the chosen criterion
% if (~isfield(options, 'flag_fienup'))
%     options.flag_fienup = false;
% end
% if (~isfield(options, 'w'))
%     options.w = ones(npix_H,npix_W);
% end
% 
% if (strcmp(options.type_obj,'dephasing') || strcmp(options.type_obj,'absorbing'))
%     crit = @(x) (critWLSlinear(x,y,Hz,H_z,-1,options.w));
% elseif (strcmp(options.type_obj,'unknown') && options.flag_fienup)
%     crit = @(x) (critFienup(x,y,Hz,H_z,-1,options.w));
% elseif (strcmp(options.type_obj,'unknown') && ~options.flag_fienup) 
%     crit = @(x) (critWLS(x,y,Hz,H_z,-1,options.w));
% end

%% Extract soft-thresholding features
if (options.mu>0.0)
    flag_softthreshod = true;
    % determine behaviour of the soft-thresholding if constraints are set
    % and TYPE_OBJ is "dephasing" or "absorbing"
    if (strcmp(options.type_obj,'dephasing'))
        if (options.iconst(1)>=0.0 && options.iconst(2)>0.0 && options.iconst(1)<=options.iconst(2))
            disp('Positivity constraint set.');
            flag_const = 1;                                                         % positivity
        elseif (options.iconst(2)<=0.0 && options.iconst(1)<0.0 && options.iconst(2)<=options.iconst(1))
            disp('Negativity constraint set.');
            flag_const = -1;                                                        % negativity
        else
            disp('No constraint set.');
            flag_const = 0;                                                         % no constraint
        end
    else
        if (options.rconst(1)>=0.0 && options.rconst(2)>0.0 && options.rconst(1)<=options.rconst(2))
            disp('Positivity constraint set.');
            flag_const = 1;                                                         % positivity
        elseif (options.rconst(2)<=0.0 && options.rconst(1)<0.0 && options.rconst(2)<=options.rconst(1))
            disp('Negativity constraint set.');
            flag_const = -1;                                                        % negativity
        else
            disp('No constraint set.');
            flag_const = 0;                                                         % no constraint
        end
    end  
else
    flag_softthreshod = false;
end

if (~isfield(options, 'support'))
    flag_support = false;
else
    if (any(size(options.support)~=[npix_W, npix_H]))
        error('SUPPORT mask must be same size as X0 and Y.');
    end
    flag_support = true;
end

if (~isfield(options, 'maxiter') || ~isscalar(options.maxiter))
    options.maxiter = 100;
end

if (~isfield(options, 'flag_cost'))
    options.flag_cost = false;
end

if (~isfield(options, 'xtrue'))
    flag_snr = false;
else
    if (any(size(options.xtrue)~=[npix_W, npix_H]))
        error('GROUND TRUTH image must be same size as X0 and Y.');
    end
    flag_snr = true;
    normxtrue = norm(options.xtrue);
    evolsnr = 20*log10(normxtrue/norm(options.xtrue(:)-x0(:)));
end

if (~isfield(options, 'flag_evolx'))
    options.flag_evolx = false;
end
if (options.flag_evolx)
    evolx = norm(x0(:));
end

if (~isfield(options, 'verbose'))
    options.verbose = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%% RECONSTRUCTION WITH FISTA ALGORITHM %%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xprev = x0;

%% INITIALIZE FISTA
uopt = xprev;
sinterp_prev = 1.0;

%% Compute first cost (data-fidelity)
if (options.flag_cost)
    [fxopt,Gxopt,c] = crit(xprev);
    evolcost = fxopt + options.mu*sum(abs(xprev(:))) ; 
end

if options.verbose
    fprintf('Iter:\t%03d\t| ', 0);
    if (options.flag_cost)
        fprintf('Cost:\t%5.2e\t| ', evolcost);
    end
    if (flag_snr)
        fprintf('SNR:\t%5.2e dB\t| ', evolsnr);
    end
    if (options.flag_evolx)
        fprintf('Evol X:\t%5.2e\t| ', evolx);
    end
    fprintf('\n');
end

%% GO ITERATE
for i=1:options.maxiter
    
    %% Gradient descent step
    [fxopt,Gxopt,c] = crit(uopt);
    uopt = uopt - (1.0/options.Lip)*Gxopt;
      
    %% Apply bound constraints  
    xopt = uopt;
    if (strcmp(options.type_obj,'dephasing'))
        idoutconst = find(xopt<options.iconst(1) | xopt>options.iconst(2));
        xopt(idoutconst) = 0.0;
    elseif (strcmp(options.type_obj,'absorbing'))
        idoutconst = find(xopt<options.rconst(1) | xopt>options.rconst(2));
        xopt(idoutconst) = 0.0;
    else
        xoptreal = xopt(:,:,1);
        idoutconst = find(xoptreal<options.rconst(1) | xoptreal>options.rconst(2));
        xoptreal(idoutconst) = 0.0;
        
        xoptimag = xopt(:,:,2);
        idoutconst = find(xoptimag<options.iconst(1) | xoptimag>options.iconst(2));
        xoptimag(idoutconst) = 0.0;
        
        xopt(:,:,1) = xoptreal;
        xopt(:,:,2) = xoptimag;
        
        clear xoptreal xoptimag;
    end
    
    %% Apply soft-thresholding operator
    if (flag_softthreshod)
        if (strcmp(options.type_obj,'dephasing') || strcmp(options.type_obj,'absorbing'))
            xopt = softThresholdingOperator(xopt,options.mu,flag_const);
        else
            xopt = softThresholdingOperator(xopt,options.mu,flag_const,true);
        end
    end
    
    %% Apply support constraint
    if (flag_support)
        xopt(:,:,1) = support.*xopt(:,:,1);
        xopt(:,:,2) = support.*xopt(:,:,2);
    end
    
    %% Interpolation step
    sinterp = 0.5*(1+sqrt(1+4*(sinterp_prev^2)));
    uopt = xopt + ((sinterp_prev-1.0)/sinterp)*(xopt - xprev);
   
    %% Compute analysis metrics
    
    % Compute Cost (data-fidelity)
    if (options.flag_cost)      
        [fxopt,Gxopt,c] = crit(xopt);
        evolcost(i+1) = fxopt + options.mu*sum(abs(xopt(:))) ;
    end
    
    % Compute SNR if required
    if (flag_snr)
        evolsnr(i+1) = 20*log10(normxtrue/norm(xtrue(:)-xopt(:)));
    end
    
    % Compute Evol X if required
    if (options.flag_evolx)
        evolx(i+1) = norm(xopt(:)-xprev(:));
    end
    
    %% Save previous iterate
    xprev = xopt;
    sinterp_prev = sinterp;
     
    %% Display convergence information
    if options.verbose
        fprintf('Iter:\t%03d\t| ', i);
        if (options.flag_cost)
            fprintf('Cost:\t%5.2e\t| ', evolcost(i+1));
        end
        if (flag_snr)
            fprintf('SNR:\t%5.2e dB\t| ', evolsnr(i+1));
        end
        if (options.flag_evolx)
            fprintf('Evol X:\t%5.2e\t| ', evolx(i+1));
        end
        fprintf('\n');
    end
    
end

if (options.flag_cost)
    varargout{1} = evolcost;
    if (flag_snr)
        varargout{2} = evolsnr;
        if (flag_evolx)
            varargout{3} = evolx;
        end
    else
        if (options.flag_evolx)
            varargout{2} = evolx;
        end
    end
else
    if (flag_snr)
        varargout{1} = evolsnr;
        if (flag_evolx)
            varargout{2} = evolx;
        end
    else
        if (flag_evolx)
            varargout{1} = evolx;
        end
    end
end

end

        
        



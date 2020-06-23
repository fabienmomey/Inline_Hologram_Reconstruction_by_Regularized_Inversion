function [xopt,evolcost,varargout] = algoFienupER(x0,crit,options)
% [xopt,varargout] = algoFienupER(x0,crit,options)
%
%   This function performs the Fienup's Error-Reduction algorithm [1]
%   for phase retrieval from intensity measurements from a transmittance
%   plane zA propagated to a plane zB.
%
%   X0: initial image guess (complex deviation from the unit transmittance
%       plane).
%
%   CRIT: function handle to the chosen criterion. It embeds:
%       Y: data image (square root of intensity measurements).
%       Hz: function handle to perform the propagation operator  (see
%           getFresnelPropagation and propagationOperator functions).
%       H_z: function handle to perform the backpropagation operator (see
%            getFresnelPropagation and propagationOperator functions).
%
%   OPTIONS: (structure of algorithm parameters)
%       * BETA       : convergence parameter for the input-output strategy
%                      (default: 1)
%                      If BETA=1, the Error-Reduction algorithm
%                       (~ Gerchberg-Saxton [2]) is performed.
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
%       * MU         :  if set, a soft-tresholding (prox L1) of parameter
%                       MU will be performed.
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
% - [1] J. Fienup, “Phase retrieval algorithms: a comparison,”
%                   Applied optics, vol. 21, no. 15, pp. 2758–2769, 1982.
% - [2] R. Gerchberg and W. Saxton, “A practical algorithm for the
%                   determination of phase from image and diffraction
%                   plane pictures,” Optik, vol. 35, p. 237, 1972.
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

%% BETA convergence relaxation parameter
%% -> if BETA = 1 (default) => Error-Reduction algorithm
%% -> if BETA < 1 => Hybrid Input-Output (HIO) algorithm (see ref [1])
if (~isfield(options, 'beta'))
    options.beta = 1.0;
end
beta = options.beta;
if (beta==1)
    warning('BETA=1: Error-Reduction algorithm is performed.');
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
%% %%%%%%%%%%%%%%%% RECONSTRUCTION WITH FIENUP ALGORITHM %%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xopt = x0;
xprev = x0;

%% Step 1 : propagation step
%xopt_prime = Hz*xopt+1.0;

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

%% Initialize the scaling factor
c=1.0;

for i=1:options.maxiter
    
    %% Gradient descent step
    [fxopt,Gxopt,c] = crit(xopt);
    xopt = xopt - 0.5*Gxopt;
    
    %% Apply bound constraints
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
        



function [xopt, fxopt, Gxopt, c, varargout] = algoRI(x0,crit,options,varargin)
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
%       Hz: function handle to perform the propagation operator (see 
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
%       * FFLAG_LINEARIZE:  a flag (default: false) for applying a 
%                           linearization of the intensity formation model
%                           (only useful if TYPE_OBJ ~= 'unknown').
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2020, Fabien Momey
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract size (in pixels) of the field of view
[npix_W, npix_H, ncomp] = size(x0);

%% Extract reconstruction parameters

%% Is a smooth regularizer given ?
flag_reg = false ;
if (isfield(options, 'reg'))
    reg = options.reg ;
    flag_reg = true ;
end

%% Lipschitz constant L
flag_backtrack = true;
if (isfield(options, 'Lip') && ~isempty(options.Lip))
    disp('FISTA without backtracking ...');
    flag_backtrack = false;
    Lip = options.Lip;
else
    disp('FISTA with backtracking ...');
    Lip = 1.0;
    eta = 1.1;
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
    % options.rconst = [max(options.rconst(1),-2.0),min(options.rconst(2),0.0)];
end

if (~isfield(options, 'iconst'))
    warning('ICONST not set: set to [-Inf,Inf] (no constraint).');
    options.iconst = [-1.0,1.0];
else
    if (any(size(options.iconst)~=[1,2]) && any(size(options.iconst)~=[2,1]))
        error('ICONST must be a 2-scalar vector.');
    end
    options.iconst=sort(options.iconst);
    % options.iconst = [max(options.iconst(1),-1.0),min(options.iconst(2),1.0)];
end

%% Extract TYPE_OBJ and adapt hard constraints bounds 
if (~isfield(options, 'type_obj'))
    options.type_obj = 'unknown';
end

if (~isfield(options, 'flag_linearize'))
    options.flag_linearize = false;
end

if (strcmp(options.type_obj,'dephasing'))
    options.rconst = [0,0];
elseif (strcmp(options.type_obj,'absorbing'))
    options.iconst = [0.0,0.0];
    options.rconst = [options.rconst(1),options.rconst(2)];
end

%% Extract soft-thresholding features
if (options.mu>0.0)
    flag_softthreshod = true;
    % determine behaviour of the soft-thresholding if constraints are set
    % and TYPE_OBJ is "dephasing" or "absorbing" and FLAG_LINEARIZE is true
    if (strcmp(options.type_obj,'dephasing') && options.flag_linearize)
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
    elseif (strcmp(options.type_obj,'absorbing') && options.flag_linearize)
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
    else
        flag_const = 0;  
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
[fuopt,Guopt,c] = crit(uopt);
if (flag_reg)
    [fureg,Gureg] = reg(uopt) ;
    fuopt = fuopt + fureg ;
    Guopt = Guopt + Gureg ;
end
fuprev = fuopt ;
Guprev = Guopt ;

%% Compute first cost (data-fidelity)
if (options.flag_cost)
    evolcost = fuopt + options.mu*sum(abs(xprev(:))) ; 
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
    ik = 0 ;
    flag_backtrack_continue = flag_backtrack ;
    while (ik<=0 || flag_backtrack_continue==true)
        %% Gradient descent step
        % [fxopt,Gxopt,c] = crit(uopt);
        % if (flag_reg)
        %     [fxreg,Gxreg] = reg(uopt) ;
        %     fxopt = fxopt + fxreg ;
        %     Gxopt = Gxopt + Gxreg ;
        % end
        uopt_new = uopt - (1.0/Lip)*Guprev;

        %% Apply bound constraints
        xopt = uopt_new;
        if (strcmp(options.type_obj,'dephasing') && options.flag_linearize)
            idoutconst = find(xopt<options.iconst(1) | xopt>options.iconst(2));
            xopt(idoutconst) = 0.0;
        elseif (strcmp(options.type_obj,'absorbing') && options.flag_linearize)
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

        %% Apply support constraint
        if (flag_support)
            xopt(:,:,1) = support.*xopt(:,:,1);
            xopt(:,:,2) = support.*xopt(:,:,2);
        end

        %% Apply soft-thresholding operator
        cost_uopt = 0.0 ;
        if (flag_softthreshod)
            if (strcmp(options.type_obj,'dephasing') || strcmp(options.type_obj,'absorbing'))
                [xopt,cost_uopt] = softThresholdingOperator(xopt,options.mu,flag_const);
            else
                [xopt,cost_uopt] = softThresholdingOperator(xopt,options.mu,flag_const,true);
            end
        end
        
        %% Re-calculate the criterion
        [fuopt,Guopt,c] = crit(xopt);
        if (flag_reg)
            [fureg,Gureg] = reg(xopt) ;
            fuopt = fuopt + fureg ;
            Guopt = Guopt + Gureg ;
        end

        %%
        d_xopt_uopt = xopt-uopt ;
        quopt = fuprev + sum(d_xopt_uopt.*Guprev,'all') + (Lip/2.0)*sum(abs(d_xopt_uopt(:)).^2) + cost_uopt ;
        %% Display backtracking information (for debugging)
        % if options.verbose
        %     fprintf('Eval. backtracking:\t%03d\t| ', ik);
        %     fprintf('fuopt:\t%5.2e\t| ', fuopt);
        %     fprintf('quopt:\t%5.2e\t| ', quopt);
        %     fprintf('L:\t%5.2e\t| ', Lip);
        %     fprintf('eta:\t%5.2e\t| ', eta^ik);
        %     fprintf('\n');
        % end
        %% Evaluate the stopping criterion of the backtracking procedure
        ik = ik+1 ;
        if (fuopt <= quopt || (~flag_backtrack_continue && ik>0))
            %% Stop backtracking flag
            flag_backtrack_continue = false ;
            uopt = uopt_new ; 
        else
            %% Update the Lipshitz constant estimate
            Lip = (eta^ik)*Lip ;
        end
    end
    
    %% Interpolation step
    sinterp = 0.5*(1+sqrt(1+4*(sinterp_prev^2)));
    uopt = xopt + ((sinterp_prev-1.0)/sinterp)*(xopt - xprev);

    %% Re-calculate the criterion
    [fuopt,Guopt,c] = crit(uopt);
    if (flag_reg)
        [fureg,Gureg] = reg(uopt) ;
        fuopt = fuopt + fureg ;
        Guopt = Guopt + Gureg ;
    end
    %% Compute analysis metrics
    
    % Compute Cost (data-fidelity)
    if (options.flag_cost)
        [fxopt,Gxopt,c] = crit(xopt);
        if (flag_reg)
            [fxreg,Gxreg] = reg(xopt) ;
            fxopt = fxopt + fxreg ;
        end
        evolcost(i+1) = fxopt + cost_uopt ;
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
    fuprev = fuopt ;
    Guprev = Guopt ;
    xprev = xopt;
    sinterp_prev = sinterp;
     
    %% Display convergence information
    if options.verbose
        fprintf('Iter:\t%03d\t| Eval. backtracking:\t%03d\t| ', i, ik);
        fprintf('L:\t%5.2e\t| ', Lip);
        fprintf('Step (1/L):\t%5.2e\t| ', 1.0/Lip);
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

%% Returned values
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

        
        



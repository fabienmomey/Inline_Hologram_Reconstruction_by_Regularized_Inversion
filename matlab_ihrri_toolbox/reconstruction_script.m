%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline Hologram Reconstruction by Regularized Inversion (IHRRI) toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This toolbox implements Inverse Problems based algorithms dedicated to
% image reconstruction in digital in-line holographic microscopy (DHM).
% Theoretical aspects on DHM and inverse approaches are developed in a
% tutorial published in JOSA A [1], and this code constitutes a 
% demonstrator of the algorithms presented in this publication.
%
% Here is the main reconstruction script that can be executed as is, and
% can perform a reconstruction from in-line hologram data whose parameters
% (data and results saving paths, calibration, algorithm settings) have to
% be set in the file "parameters.m" (refer to it for more details).
%
% All settings are stored in a global structure EXPE which also store
% the reconstruction results. At the end of the script, this structure is
% saved in a MAT-file "experiment.m" in a timstamped subdirectory contained
% in the "results" directory that is named as the data hologram file.
%
% Summary
% -------
%
% The code is able to perform 2 kind of "Inverse problems" algorithms
% aiming at reconstructing an image X from an intensity in-line hologram 
% image Y. In this code, X is 2-component image ([width,height,2]), each
% one corresponding respectively to the real and imaginary part of the
% complex deviation from the unit transmittance plane :
%
%                   T = 1 + X(:,:,1) + i X(:,:,2)
%
% The 2 implemented reconstruction algorithm are:
%   - The Fienup's Error-Reduction algorithm [2], which stands for a
%   gradient descent algorithm aiming at solving the following
%   problem:
%   
%   X   = ARG MIN   || C*M(X) - Y ||_2^2        s.t.    X in Omega
%            x
%
%   where 
%       - M(X) = |1 + H.X| is the square root intensity hologram formation 
%       model of X using a convolutive propagation operator H (the Fresnel 
%       kernel is available in this toolbox).
%       - Omega stands for the validity domain of X that is enforced as a
%       "projection on constraints" operator in the algorithm. The domain
%       Omega takes the form of bound constraints applied on the real and
%       imaginary parts of X:
%           * RCONST = [XRMIN,XRMAX]
%           * ICONST = [XIMIN,XIMAX]
%
%       - C is a scaling factor that accounts for the intensity of the 
%       incident wave |a_0|^2 as well as the detector gain and quantum 
%       efficiency [1].
%
%   - The FISTA algorithm [3], which is a proximal gradient descent
%   algorithm aiming at solving the following sparsity problem:
%   
%   X   = ARG MIN   || C*M(X) - Y ||_W^2 + mu * || X ||_1
%            x
%
%                                          s.t.    X in Omega
%
%   where 
%       - M(X) is still the intensity hologram formation model which can 
%       take 2 forms:
%           * Considering purely and weakly dephasing or purely absorbing 
%           objects, we can use a linearized intensity hologram formation
%           model:
%
%           M(X) = 1 + G.X ~ |1 + H.X|^2 
%
%           where X becomes purely real image (size [width,height]) and G
%           is a purely real kernel (see the code's documentation for
%           details) which depends on the convolutive propagation operator 
%           H.
%
%           * Considering a unknown object:
%
%           M(X) = |1 + H.X|^2 
%           
%           that models the "full" intensity hologram image from the
%           "complex" image X.
%
%       - Omega stands for the validity domain of X that is enforced as a
%       "projection on constraints" operator in the algorithm.
%       - C is a scaling factor that accounts for the intensity of the 
%       incident wave |a_0|^2 as well as the detector gain and quantum 
%       efficiency [1].
%       - W is the inverse noise covariance matrix C^{-1}.
%
%       - || X ||_1 is a regularizer enforcing a sparsity constraint,
%       weighted by a parameter MU, and applied as a soft-thresholding
%       operator in the algorithm.
%
% Note that the code is able to perform fied-of-view extension in a
% rigorous way by setting an extension factor in the settings (see
% "parameters.m").
%
% In the future, an edge-preserving regularizer will be implemented which
% can be added to the reconstruction criterion (see [1] for theoretical
% details).
%
% References
%
% - [1] F. Momey, L. Denis, T. Olivier, C. Fournier, "From Fienup’s phase 
%                   retrieval techniques to regularized inversion for 
%                   in-line holography: tutorial," JOSA A, vol. 36, no. 12, 
%                   D62-D80, 2019. 
% - [2] J. Fienup, “Phase retrieval algorithms: a comparison,”
%                   Applied optics, vol. 21, no. 15, pp. 2758–2769, 1982.
% - [3]  A. Beck and M. Teboulle, “Fast gradient-based algorithms for con-
%       strained total variation image denoising and deblurring problems,”
%       IEEE Trans. Image Process. 18, 2419–2434, 2009.
%
% Created: 04/06/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516, 
%           Université Jean Monnet, 
%           F-42000 Saint-Étienne, 
%           France
%           fabien.momey@univ-st-etienne.fr
%
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

% clear all ;
close all ;
clc ;

addpath(genpath('./'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD EXPERIMENT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%% LOADING DATA %%
%%%%%%%%%%%%%%%%%%
data = double(imread([EXPE.holodir_data,EXPE.holodatafile]));
data = data/median(data(:));

% Get hologram dimensions
EXPE.fov_width = size(data,2) ;           % field-of-view width in pixels
EXPE.fov_height = size(data,1) ;          % field-of-view height in pixels

if (EXPE.flag_display)
    ihrri_show(data, 'Data');
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.data = data;
clear data;

%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL CALIBRATION %%
%%%%%%%%%%%%%%%%%%%%%%%

EXPE.flag_fov_extension = false;
if (EXPE.fov_extension_factor>1.0)
    EXPE.flag_fov_extension = true;
    % Adapt the fov size with hte extension factor
    newfovsize = size(fovExtensionOperator(zeros(EXPE.fov_height,EXPE.fov_width),EXPE.fov_extension_factor));
    EXPE.fov_width = newfovsize(2);
    EXPE.fov_height = newfovsize(1);
    clear newfovsize;
end

%% Get Propagator Kernel
if (strcmp(EXPE.type_prop,'Fr')) % Use Fresnel propagator
    if EXPE.flag_pad
        [Hz]=getFresnelPropagator(2*EXPE.fov_width,...
            2*EXPE.fov_height,...
            EXPE.pixel_size,...
            EXPE.lambda,...
            EXPE.n_0,...
            EXPE.z_s,...
            EXPE.NA);
        if (~strcmp(EXPE.type_obj,'unknown') && EXPE.flag_linearize)
            [Gz]=getFresnelPropagator(2*EXPE.fov_width,...
                2*EXPE.fov_height,...
                EXPE.pixel_size,...
                EXPE.lambda,...
                EXPE.n_0,...
                EXPE.z_s,...
                EXPE.NA,...
                false,...
                'intensity',... % this code only deals with intensity holograms
                EXPE.flag_linearize,...
                EXPE.type_obj);
            EXPE.Gz = Gz;
        end
    else
        [Hz]=getFresnelPropagator(EXPE.fov_width,...
            EXPE.fov_height,...
            EXPE.pixel_size,...
            EXPE.lambda,...
            EXPE.n_0,...
            EXPE.z_s,...
            EXPE.NA);
        if (~strcmp(EXPE.type_obj,'unknown') && EXPE.flag_linearize)
            [Gz]=getFresnelPropagator(EXPE.fov_width,...
                EXPE.fov_height,...
                EXPE.pixel_size,...
                EXPE.lambda,...
                EXPE.n_0,...
                EXPE.z_s,...
                EXPE.NA,...
                false,...
                'intensity',... % this code only deals with intensity holograms
                EXPE.flag_linearize,...
                EXPE.type_obj);
            EXPE.Gz = Gz;
        end
    end
elseif(strcmp(EXPE.type_prop,'RS')) % Use Rayleigh-Sommerfeld propagator
    if EXPE.flag_pad
        [Hz]=getRayleighSommerfeldPropagator(2*EXPE.fov_width,...
            2*EXPE.fov_height,...
            EXPE.pixel_size,...
            EXPE.lambda,...
            EXPE.n_0,...
            EXPE.z_s, ...
            EXPE.NA);
        if (~strcmp(EXPE.type_obj,'unknown') && EXPE.flag_linearize)
            [Gz]=getRayleighSommerfeldPropagator(2*EXPE.fov_width,...
                2*EXPE.fov_height,...
                EXPE.pixel_size,...
                EXPE.lambda,...
                EXPE.n_0,...
                EXPE.z_s,...
                EXPE.NA,...
                false,...
                'intensity',... % this code only deals with intensity holograms
                EXPE.flag_linearize,...
                EXPE.type_obj);
            EXPE.Gz = Gz;
        end
    else
        [Hz]=getRayleighSommerfeldPropagator(EXPE.fov_width,...
            EXPE.fov_height,...
            EXPE.pixel_size,...
            EXPE.lambda,...
            EXPE.n_0,...
            EXPE.z_s, ...
            EXPE.NA);
        if (~strcmp(EXPE.type_obj,'unknown') && EXPE.flag_linearize)
            [Gz]=getRayleighSommerfeldPropagator(EXPE.fov_width,...
                EXPE.fov_height,...
                EXPE.pixel_size,...
                EXPE.lambda,...
                EXPE.n_0,...
                EXPE.z_s,...
                EXPE.NA,...
                false,...
                'intensity',... % this code only deals with intensity holograms
                EXPE.flag_linearize,...
                EXPE.type_obj);
            EXPE.Gz = Gz;
        end
    end
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.Hz = Hz;
clear Gz Hz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMPLE BACK-PROPAGATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xbackprop = propagationOperator(sqrt(EXPE.data)-1.0,conj(EXPE.Hz));
% 
% if (EXPE.flag_display)
%     ihrri_show(abs(xbackprop), 'Backpropagation (modulus)');
%     ihrri_show(angle(xbackprop), 'Backpropagation (phase)');
% end
% 
% % FILL EXPERIMENT REPORT STRUCT
% EXPE.xbackprop = xbackprop;
% clear xbackprop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERIZING RECONSTRUCTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZATION PARAMETERS  
if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
    && EXPE.flag_linearize)
    EXPE.o0 = zeros([EXPE.fov_height, EXPE.fov_width]);
else
    EXPE.o0 = zeros([EXPE.fov_height, EXPE.fov_width, 2]);
end

% GRADIENT IMAGE FUNCTION HANDLES (FOR EDGE-PRESERVING SMOOTHING REGULARIZATION)
if (isfield(EXPE,'muEdgePres') && any(EXPE.muEdgePres ~= 0.0))
    if EXPE.flag_pad
        EXPE.Grad = getGradientOperator(2*EXPE.fov_width, 2*EXPE.fov_height) ;
    else
        EXPE.Grad = getGradientOperator(EXPE.fov_width, EXPE.fov_height) ;
    end
    % Gradient in X direction
    Gx = @(o) (propagationOperator(o,EXPE.Grad,EXPE.flag_pad)) ;
    % Gradient in Y direction
    Gy = @(o) (propagationOperator(o,conj(EXPE.Grad'),EXPE.flag_pad)) ;
    % Adjoint gradient in X direction
    GxT = @(o) (propagationOperator(o,conj(EXPE.Grad),EXPE.flag_pad)) ;
    % Adjoint gradient in Y direction
    GyT = @(o) (propagationOperator(o,EXPE.Grad',EXPE.flag_pad)) ;
    % If epsilonEdgePres is not given, it is arbitrarily set to a very smal number
    if (~isfield(EXPE,'epsilonEdgePres'))
        EXPE.epsilonEdgePres = 1.0e-12 ;
    end
end

% PROPAGATION AND BACKPROPAGATION FUNCTION HANDLES
if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
        && EXPE.flag_linearize)
    if (~EXPE.flag_fov_extension)
        Propag = @(o) (propagationOperator(o,EXPE.Gz,EXPE.flag_pad));
        BackPropag = @(o) (propagationOperator(o,EXPE.Gz,EXPE.flag_pad));
    else
        Propag = @(o) (fovExtensionOperator(...
            propagationOperator(o,EXPE.Gz,EXPE.flag_pad),...
            EXPE.fov_extension_factor,...
            true));
        BackPropag = @(o) (propagationOperator(...
            fovExtensionOperator(o,EXPE.fov_extension_factor),...
            EXPE.Gz,...
            EXPE.flag_pad));
    end
else
    if (~EXPE.flag_fov_extension)
        Propag = @(o) (propagationOperator(o,EXPE.Hz,EXPE.flag_pad));
        BackPropag = @(o) (propagationOperator(o,conj(EXPE.Hz),EXPE.flag_pad));
    else
        Propag = @(o) (fovExtensionOperator(...
            propagationOperator(o,EXPE.Hz,EXPE.flag_pad),...
            EXPE.fov_extension_factor,...
            true));
        BackPropag = @(o) (propagationOperator(...
            fovExtensionOperator(o,EXPE.fov_extension_factor),...
            conj(EXPE.Hz),...
            EXPE.flag_pad));
    end
end

% GET THE CRITERION
if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
        && EXPE.flag_linearize)
    Crit = @(x) (critWLSlinear(x,EXPE.data,Propag,BackPropag,-1));  
elseif ((strcmp(EXPE.type_obj,'unknown') && EXPE.flag_fienup) || strcmp(EXPE.flag_rec_meth,'Fienup'))
    Crit = @(x) (critFienup(x,sqrt(EXPE.data),Propag,BackPropag,-1));
else
    Crit = @(x) (critWLS(x,EXPE.data,Propag,BackPropag,-1)) ;  
end
% EDGE-PRESERVING REGULARIZATION
if (isfield(EXPE,'muEdgePres') && any(EXPE.muEdgePres ~= 0.0))
    if (strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))
        Reg = @(x) (critEdgePreservingSmoothing(x,Gx,Gy,GxT,GyT,EXPE.muEdgePres,EXPE.epsilonEdgePres)) ;
    else
        Reg = @(x) (critEdgePreservingSmoothing(x,Gx,Gy,GxT,GyT,EXPE.muEdgePres,EXPE.epsilonEdgePres,true)) ;
    end
end


% PREPARE ALGORITHM
switch EXPE.flag_rec_meth
    case 'RI'
        % Get Lipschitz constant if backtracking is not set on
        if (~EXPE.flag_backtrack)
            if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
                    && EXPE.flag_linearize)
                EXPE.Lip = 2*max(EXPE.Gz(:).*EXPE.Gz(:));
            else
                EXPE.Lip = 2*max(conj(EXPE.Hz(:)).*EXPE.Hz(:))^2;
            end
        else
            EXPE.Lip = [] ;
        end
          
        RECoptions = struct('type_obj',EXPE.type_obj,...
            'flag_linearize',EXPE.flag_linearize,...
            'rconst',EXPE.real_constraint,...
            'iconst',EXPE.imag_constraint,...
            'mu',EXPE.muSparse,...
            'maxiter', EXPE.maxiter,...
            'flag_cost',true,...
            'flag_evolx',false,...
            'verbose',true,...
            'Lip',EXPE.Lip);
    case 'Fienup'
        RECoptions = struct('type_obj',EXPE.type_obj,...
            'rconst',EXPE.real_constraint,...
            'iconst',EXPE.imag_constraint,...
            'mu',EXPE.muSparse,...
            'maxiter', EXPE.maxiter,...
            'flag_cost',true,...
            'flag_evolx',false,...
            'verbose',true);
end
if (isfield(EXPE,'muEdgePres') && any(EXPE.muEdgePres ~= 0.0))
    RECoptions.('reg') = Reg ;
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.RECoptions = RECoptions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAUNCH RECONSTRUCTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Let''s reconstruct !');
switch EXPE.flag_rec_meth
    case 'RI'
        [RECxopt,RECfxopt,RECGxopt,RECc,RECevolcost] = algoRI(EXPE.o0,Crit,RECoptions);
    case 'Fienup'
        [RECxopt,RECfxopt,RECGxopt,RECc,RECxopt,RECevolcost] = algoFienupER(EXPE.o0,Crit,RECoptions);
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.evolcost = RECevolcost;
EXPE.xopt = RECxopt;
EXPE.fxopt = RECfxopt;
EXPE.Gxopt = RECGxopt;
EXPE.c = RECc;

% DISPLAY RECONSTRUCTION
if (EXPE.flag_display)
    %% Reconstruction
    if (~strcmp(EXPE.type_obj,'Fienup') && ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
                && EXPE.flag_linearize))
        if (strcmp(EXPE.type_obj,'dephasing'))
            ihrri_show(RECxopt,'Reconstructed phase');
        elseif (strcmp(EXPE.type_obj,'absorbing'))
            ihrri_show(-RECxopt,'Reconstructed opacity');
        end
    else
        RECxopt = 1.0 + RECxopt(:,:,1) + 1i * RECxopt(:,:,2);
        ihrri_show(angle(RECxopt),'Reconstructed phase');
        ihrri_show(abs(RECxopt),'Reconstructed modulus');
    end
    %% Residues
    [fxopt,gxopt,c,residues] = Crit(EXPE.xopt);
    ihrri_show(residues,'Residues');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS UNSTACKING AND SAVING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVING
save([EXPE.holodir_results_timestamp,'experiment.mat'],'-struct','EXPE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline Hologram Reconstruction by Regularized Inversion (IHRRI) toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script allows to set all the parameters (data and results saving 
% paths, calibration, algorithm settings) to perform a reconstruction of a
% given in-line hologram.
%
% All settings are stored in a global structure EXPE which will also save
% all the reconstruction results.
%
% Here is a list of all the parameters:
% - HOLODIR:        working directory (default: the current directory where
%                   this code is installed).
% - HOLODIR_DATA:   data directory (the DATA directory has to be located
%                   at the root of HOLODIR directory. Some test holograms
%                   are already available in this directory.
% - HOLODATAFILE:   name of the data file.
% - HOLODIR_RESULTS:results directory that will be created at the root of 
%                   HOLODIR directory. Then timestamped subdirectories for 
%                   the current experiment are automatically created in the 
%                   folder tree from HOLODIR_RESULTS.
%
% - FLAG_REC_METH:  'Fienup' or 'RI' => Reconstruction method.
% - FLAG_FIENUP:    true or false => Choose Fienup criterion 
%                   (only useful if flag_rec_meth = 'RI').
% - TYPE_OBJ:       'unknown', 'dephasing' or 'absorbing' => the object of 
%                   interest is purely 'dephasing' or 'absorbing', 
%                   or 'unknown'. It allows to define the propagation
%                   kernel and default bound constraints.
% - FLAG_LINEARIZE: true or false => linearization of the intensity 
%                   formation model (only useful if TYPE_OBJ ~= 'unknown').
%                   If flag_rec_meth = 'Fienup', this flag must be set to 
%                   false to get the correct propagation kernel.
% - FLAG_PAD:       true (recommended) or false => allows zero-padding for  
%                   performing convolutions.
% - FLAG_DISPLAY:   true or false => allows to display results.
%
% - Z_S:            Distance (in m) from the sensor plane to the object 
%                   plane.
% - MAG:            Lens magnification.
% - LAMBDA:         Illumination wavelength (in m).
% - N_0:            Medium refractive index.
% - PIXEL_SIZE:     pixel size (both in x and y) (in m).
% - FOV_WIDTH:      field-of-view width (in pixels).
% - FOV_HEIGHT:     field-of-view height (in pixels).
% - FOV_EXTENSION_FACTOR:   field-of-view extension factor (cannot be <1 ; 
%                           if 1 => no fov extension).
% - REAL_CONSTRAINT:a 2-element vector giving hard constraint parameter for 
%                   the real part of X
%                           \_ TYPE_OBJ = 'dephasing'
%                               \_ default: [0,0] (X is purely imaginary)
%                           \_ TYPE_OBJ = 'absorbing'
%                               \_ default: [-1,0] (X is purely real)
%                           \_ TYPE_OBJ = 'unknown'
%                               \_ default: [-2,0] (because 0 < |T| < 1
%                                                   and -1 < cos(phi) < 1)
% - IMAG_CONSTRAINT:a 2-element vector giving hard constraint
%                          parameter for the imag part of X
%                           \_ TYPE_OBJ = 'dephasing'
%                               \_ default: [-1,1] (X is purely imaginary)
%                           \_ TYPE_OBJ = 'absorbing'
%                               \_ default: [0,0] (X is purely real)
%                           \_ TYPE_OBJ = 'unknown'
%                               \_ default: [-1,1] (because 0 < |T| < 1
%                                                   and -1 < sin(phi) < 1)
% - MUSPARSE:       hyperparameter for the sparsity constraint 
%                   (soft-thresholding operator)
% - MUEDGEPRES:     hyperparameter \mu for the edge-preserving regularizer. 
%                   Not yet implemented.
% - EPSILONEDGEPRES:hyperparameter \epsilon for the edge-preserving regularizer. 
%                   Not yet implemented.
% - MAXITER:        Maximum number of iterations.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE EXPERIMENT REPORT STRUCT AND FILL FIRST PARAMETERS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPE = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIRECTORIES AND FILENAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global directory
EXPE.holodir = ['./'];
% Data directory
EXPE.holodir_data = [EXPE.holodir,'data/2019_06_07_billes1mu_63_JOSAA/'];
% Data filename
EXPE.holodatafile = 'Basler daA1920-30um (22030948)_20190607_122251706_0012_crop.tiff';

% Creation of the global results' directory
holodir_results = [EXPE.holodir,'results/'];
if (~exist(holodir_results,'dir'))
    mkdir(EXPE.holodir,'results/');
end

% Creation of the experiment results' directory and filenames
idext = strfind(EXPE.holodatafile,'.tiff');
EXPE.holoresultfile = EXPE.holodatafile(1:(idext-1));
EXPE.holodir_results_expe = [holodir_results,EXPE.holoresultfile,'/'];
if (~exist(EXPE.holodir_results_expe,'dir'))
    mkdir(holodir_results,EXPE.holoresultfile);
end

% Create the time stamped experiment directory
newdate_expe = datestr(now,30);
mkdir(EXPE.holodir_results_expe,newdate_expe);
EXPE.holodir_results_timestamp = [EXPE.holodir_results_expe,newdate_expe,'/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLAGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I choose the reconstruction method: 'Fienup' or 'RI'
EXPE.flag_rec_meth = 'RI';
%% Choose Fienup criterion (only useful if flag_rec_meth = 'RI')
EXPE.flag_fienup = false;

%% My object of interest is purely 'dephasing' or 'absorbing', or 'unknown'.
%% It allows to define the propagation kernel and default bound constraints.
EXPE.type_obj = 'dephasing';

%% I want a linearization of the intensity formation model: true or false
%% (only useful if type_obj != 'unknown')
%% If Fienup ER is performed, this flag must be set to false
EXPE.flag_linearize = true;

%% My reconstruction requires:
% - zero-padding for performing convolutions
EXPE.flag_pad = true;

%% I want to display results
EXPE.flag_display = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION PARAMETERS:

% - INSTRUMENTAL
EXPE.z_s = 7.2822e-06 ;          % (m) Distance from the sensor plane to 
                                 % the object plane
EXPE.mag = 56.7 ;                % Lens magnification
EXPE.lambda = 532e-9  ;          % (m) wavelength
EXPE.n_0 = 1.52 ;                % Medium refractive index (not mandatory)

% - DIGITAL
EXPE.pixel_size = 2.2e-6/EXPE.mag ;	% (m) % pixel size
EXPE.fov_width = 512 ;           % field-of-view width in pixels
EXPE.fov_height = 512 ;          % field-of-view	height in pixels
EXPE.fov_extension_factor = 1.5; % field-of-view extension factor
                            % (cannot be <1 ; if =1 => no fov extension)
% RECONSTRUCTION PARAMETERS
EXPE.real_constraint = [0.0,0.0];   % a 2-element vector giving hard constraint
%                          parameter for the real part of X
%                           \_ TYPE_OBJ = 'dephasing'
%                               \_ default: [0,0] (X is purely imaginary)
%                           \_ TYPE_OBJ = 'absorbing'
%                               \_ default: [-1,0] (X is purely real)
%                           \_ TYPE_OBJ = 'unknown'
%                               \_ default: [-2,0] (because 0 < |T| < 1
%                                                   and -1 < cos(phi) < 1)
EXPE.imag_constraint = [0.0,Inf];   % a 2-element vector giving hard constraint
%                          parameter for the imag part of X
%                           \_ TYPE_OBJ = 'dephasing'
%                               \_ default: [-1,1] (X is purely imaginary)
%                           \_ TYPE_OBJ = 'absorbing'
%                               \_ default: [0,0] (X is purely real)
%                           \_ TYPE_OBJ = 'unknown'
%                               \_ default: [-1,1] (because 0 < |T| < 1
%                                                   and -1 < sin(phi) < 1)
EXPE.muSparse = 0.01;             % hyperparameter for the sparsity constraint
                            % (soft-thresholding operator)
EXPE.muEdgePres = 0.1;           % hyperparameter \mu for the edge-preserving
                            % regularizer (if required)
                            %% NOT YET AVAILABLE
EXPE.epsilonEdgePres = 1.0e-2;   % hyperparameter \epsilon for the
                            % edge-preserving regularizer (if required) 
                            %% NOT YET AVAILABLE

% OPTIMIZATION PARAMETERS
EXPE.maxiter = 50;              % maximum number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

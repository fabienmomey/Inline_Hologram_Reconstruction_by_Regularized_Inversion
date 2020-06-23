%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters ...
%
% Licence ...
%
% Created: 04/06/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516, 
%           Université Jean Monnet, 
%           F-42000 Saint-Étienne, 
%           France
%           fabien.momey@univ-st-etienne.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE EXPERIMENT REPORT STRUCT AND FILL FIRST PARAMETERS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPE = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIRECTORIES AND FILENAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global directory
EXPE.holodir = ['/home/mof03096/Documents/Laboratoire_Hubert_Curien/Recherche/Codes_et_Donnees/',...
'Inline_Hologram_Reconstruction_by_Regularized_Inversion/matlab_ihrri_toolbox/'];
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
EXPE.holodir_results_timestamp = [EXPE.holodir_results_expe,newdate_expe];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLAGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% I choose the reconstruction method: 'Fienup' or 'RI'
EXPE.flag_rec_meth = 'Fienup';
%% Choose Fienup criterion (only useful if flag_rec_meth = 'RI')
EXPE.flag_fienup = false;

%% My object of interest is purely 'dephasing' or 'absorbing', or 'unknown':
EXPE.type_obj = 'unknown';

%% I want a linearization of the intensity formation model: true or false
%% (only useful if type_obj != 'unknown')
EXPE.flag_linearize = true;

%% My reconstruction requires:
% - zero-padding for performing convolutions
EXPE.flag_pad = true;
% % - Hard constraint: positivity (1) or negativity (-1)
% EXPE.flag_hard_constraint = 1;
% % - Edge-preserving regularization
% EXPE.flag_edge_preserving = true;

%% I want to display results
EXPE.flag_display = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION PARAMETERS:

% - INSTRUMENTAL
EXPE.z_s = 7.2822e-06 ;          % (m) Distance of the sensor of the plane
EXPE.mag = 56.7 ;                % Lens magnification
EXPE.lambda = 532e-9  ;          % (m) wavelength
EXPE.n_0 = 1.52 ;                % Medium refractive index (not mandatory)

% - DIGITAL
EXPE.pixel_size = 2.2e-6/EXPE.mag ;	% (m) % pixel size
EXPE.fov_width = 512 ;           % field-of-view width in pixels
EXPE.fov_height = 512 ;          % field-of-view	height in pixels
EXPE.fov_extension_factor = 1.0; % field-of-view extension factor
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
EXPE.muSparse = 0.0;             % hyperparameter for the sparsity constraint
                            % (soft-thresholding operator)
EXPE.muEdgePres = 0.1;           % hyperparameter \mu for the edge-preserving
                            % regularizer (if required)
EXPE.epsilonEdgePres = 1.0e-2;   % hyperparameter \epsilon for the
                            % edge-preserving regularizer (if required)       

% OPTIMIZATION PARAMETERS
EXPE.stepfista = 0.05;           % step size for FISTA algorithm
EXPE.maxiter = 10;              % maximum number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

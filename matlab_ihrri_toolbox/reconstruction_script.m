%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implements ...
%
% Licence ...
%
% References
%
% - [1] F. Momey, L. Denis, T. Olivier, C. Fournier, "From Fienup’s phase 
%                   retrieval techniques to regularized inversion for 
%                   in-line holography: tutorial," JOSA A, vol. 36, no. 12, 
%                   D62-D80, 2019. 
%
% Created: 04/06/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516, 
%           Université Jean Monnet, 
%           F-42000 Saint-Étienne, 
%           France
%           fabien.momey@univ-st-etienne.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ;
close all ;
clc ;

addpath(genpath('./'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD EXPERIMENT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if EXPE.flag_pad
    [Hz]=getFresnelPropagator(2*EXPE.fov_width,...
        2*EXPE.fov_height,...
        EXPE.pixel_size,...
        EXPE.lambda,...
        EXPE.n_0,...
        EXPE.z_s);
    if (~strcmp(EXPE.type_obj,'unknown') && EXPE.flag_linearize)
        [Gz]=getFresnelPropagator(2*EXPE.fov_width,...
            2*EXPE.fov_height,...
            EXPE.pixel_size,...
            EXPE.lambda,...
            EXPE.n_0,...
            EXPE.z_s,...
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
        EXPE.z_s);
    if (~strcmp(EXPE.type_obj,'unknown') && EXPE.flag_linearize)
        [Gz]=getFresnelPropagator(EXPE.fov_width,...
            EXPE.fov_height,...
            EXPE.pixel_size,...
            EXPE.lambda,...
            EXPE.n_0,...
            EXPE.z_s,...
            false,...
            'intensity',... % this code only deals with intensity holograms
            EXPE.flag_linearize,...
            EXPE.type_obj);
        EXPE.Gz = Gz;
    end
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.Hz = Hz;
clear Gz Hz;

%%%%%%%%%%%%%%%%%%
%% LOADING DATA %%
%%%%%%%%%%%%%%%%%%
data = double(imread([EXPE.holodir_data,EXPE.holodatafile]));
data = data/median(data(:));

% if (EXPE.flag_fov_extension)
%     data = fovExtensionOperator(data,EXPE.fov_extension_factor);
% end

if (EXPE.flag_display)
    ihrri_show(data, 'Data');
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.data = data;
clear data;

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
    EXPE.o0 = zeros([EXPE.fov_width, EXPE.fov_height]);
else
    EXPE.o0 = zeros([EXPE.fov_width, EXPE.fov_height, 2]);
end

% PROPAGATION AND BACKPROPAGATION FUNCTION HANDLES
if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
        && EXPE.flag_linearize)
    if (~EXPE.flag_fov_extension)
        Propag = @(o) (propagationOperator(o,EXPE.Gz,true));
        BackPropag = Propag;
    else
        Propag = @(o) (fovExtensionOperator(...
            propagationOperator(o,EXPE.Gz,true),...
            EXPE.fov_extension_factor,...
            true));
        BackPropag = @(o) (propagationOperator(...
            fovExtensionOperator(o,EXPE.fov_extension_factor),...
            EXPE.Gz,...
            true));
    end
else
    if (~EXPE.flag_fov_extension)
        Propag = @(o) (propagationOperator(o,EXPE.Hz,true));
        BackPropag = @(o) (propagationOperator(o,conj(EXPE.Hz),true));
    else
        Propag = @(o) (fovExtensionOperator(...
            propagationOperator(o,EXPE.Hz,true),...
            EXPE.fov_extension_factor,...
            true));
        BackPropag = @(o) (propagationOperator(...
            fovExtensionOperator(o,EXPE.fov_extension_factor),...
            conj(EXPE.Hz),...
            true));
    end
end

% GET THE CRITERION
if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
        && EXPE.flag_linearize)
    Crit = @(x) (critWLSlinear(x,EXPE.data,Propag,BackPropag,-1));
elseif ((strcmp(EXPE.type_obj,'unknown') && EXPE.flag_fienup) || strcmp(EXPE.flag_rec_meth,'Fienup'))
    Crit = @(x) (critFienup(x,sqrt(EXPE.data),Propag,BackPropag,-1));
elseif (strcmp(EXPE.type_obj,'unknown') && ~EXPE.flag_fienup)
    Crit = @(x) (critWLS(x,EXPE.data,Propag,BackPropag,-1));
end

% PREPARE ALGORITHM
switch EXPE.flag_rec_meth
    case 'RI'
        % Get Lipschitz constant
        if ((strcmp(EXPE.type_obj,'dephasing') || strcmp(EXPE.type_obj,'absorbing'))...
                && EXPE.flag_linearize)
            EXPE.Lip = 2*max(EXPE.Gz(:).*EXPE.Gz(:));
        else
            EXPE.Lip = 2*max(conj(EXPE.Hz(:)).*EXPE.Hz(:));
        end
          
        RECoptions = struct('Lip',EXPE.Lip,...
            'type_obj',EXPE.type_obj,...
            'rconst',EXPE.real_constraint,...
            'iconst',EXPE.imag_constraint,...
            'mu',EXPE.muSparse,...
            'maxiter', EXPE.maxiter,...
            'flag_cost',true,...
            'flag_evolx',false,...
            'verbose',true);
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

% FILL EXPERIMENT REPORT STRUCT
EXPE.RECoptions = RECoptions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAUNCH RECONSTRUCTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Let''s reconstruct !');
switch EXPE.flag_rec_meth
    case 'RI'
        [RECxopt,RECevolcost] = algoRI(EXPE.o0,Crit,RECoptions);
    case 'Fienup'
        [RECxopt,RECevolcost] = algoFienupER(EXPE.o0,Crit,RECoptions);
end

% FILL EXPERIMENT REPORT STRUCT
EXPE.evolcost = RECevolcost;
EXPE.xopt = RECxopt;

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

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS UNSTACKING AND SAVING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVING
switch flag_rec_meth
    case 'RI'
        holoresultfile = sprintf('%s%s_expe-%s_meth-%s_objtype-%s_pad-%d_ext-%d_oneshot-%d_Lip-%2.2e_gam-%2.2e_mu-%2.2e.mat',...
            holodir_results_expe,...
            holoresultfile,...
            flag_expe,...
            flag_rec_meth,...
            type_obj,...
            flag_pad,...
            flag_ext,...
            flag_oneshot,...
            Lip,...
            gamREC/Lip,...
            muREC);
    case 'Fienup'
        holoresultfile = sprintf('%s%s_expe-%s_meth-%s_objtype-%s_pad-%d_ext-%d_oneshot-%d_Lip-%2.2e_gam-%2.2e_beta-%2.2e_prox-%d_mu-%2.2e_support-%d.mat',...
                holodir_results_expe,...
                holoresultfile,...
                flag_expe,...
                flag_rec_meth,...
                type_obj,...
                flag_pad,...
                flag_ext,...
                flag_oneshot,...
                Lip,...
                gamREC,...
                betaREC,...
                applyproxHIO,...
                flag_support);
end
save(holoresultfile,'-struct','EXPE');





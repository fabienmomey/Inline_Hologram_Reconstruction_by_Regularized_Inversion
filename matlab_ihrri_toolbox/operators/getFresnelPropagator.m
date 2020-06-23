function [Hz] = getFresnelPropagator(npix_W, npix_H, pixel_size, lambda, n_0, z, varargin)
% [Hz] = getFresnelPropagator(npix_W, npix_H, pixel_size, lambda, n_0, z, varargin)
%
%   This function calculates the Fresnel propagation kernel used for 
%   Fresnel propagation (convolution).
%
%   The kernel Hz and its backpropagation counterpart (adjoint operator)
%   requires the parameters:
%   - NPIX_W: output kernel width in pixels.
%   - NPIX_W: output kernel height in pixels.
%   - PIXEL_SIZE: pixel size (supposed to be the same in width and height).
%   - LAMBDA: illumination wavelength.
%   - N_0: refractive index of the medium.
%   - Z: propagation distance.
% 
%   The kernel calculation is based on a particular value :
%
%       z_shannon = max([npix_W, npix_H])*(pixel_size^2)/lambda
%
%   If z > z_shannon
%   | --> Analytic impulse response in the spatial domain + FFT
%   Else
%   | --> Analytic transfert function in the Fourier domain
%   End
%
%   In VARGARIN, 3 parameters can be given:
%   - FLAG_PHASEREF: if true (default: false), the phase origin is taken 
%   into account in the calculation of the propagation kernel, 
%   that is to say that the absolute phase term exp(2i pi Z / LAMBDA) 
%   induced by the propagation distance Z is calculated in the kernel.
%
%   - TYPE_HOLOGRAM: indicates if complex wave or intensity data will be
%   used. Possible values are
%       > 'complex' (default): the Fresnel propagation kernel HZ is used in
%       the context of complex wave measurements. It requires of full
%       complex Fresnel kernel.
%
%           U = |A0|.(1 + HZ*O)
%
%           where U is the total complex wave and O represents the targeted 
%           complex deviation from the unit transmittance plane.
%
%       > 'intensity': the Fresnel propagation kernel is used in the
%       context of intensity measurements. Hence, a "adapted" Fresnel
%       kernel can be calculated depending on the TYPE_MODEL parameter.
%
%           I = |U|^2 = |A0|^2.|1 + HZ*O|^2
%
%           where I is the intensity of the total complex wave and O still 
%           represents the targeted complex deviation from the unit 
%           transmittance plane.
%
%   - FLAG_LINEARIZE: to set only if TYPE_HOLOGRAM='intensity'. It indicates 
%   which kind of model will be applied to possibly adapt the calculation 
%   of the kernel. Possible values are:
%       > false (default): the full propagation model has to be 
%       calculated. It requires the calculation of the full complex Fresnel 
%       kernel HZ.
%       > true: used in the context of TYPE_HOLOGRAM='intensity' only.
%       A linearized model of intensity measurement can be applied. It
%       yields a modified/simplified Fresnel kernel depending on the
%       TYPE_OBJ parameter. In this case, FLAG_PHASEREF will be set to
%       false.
%   - TYPE_OBJ: to set only if TYPE_LINEARIZE=true. Possible values are 
%       > 'dephasing' (default): purely and weakly dephasing object. The
%       linearized intensity measurement formation model writes [1]:
%
%           I = 1 + GZ * O
%
%           where GZ = -2 Im[HZ] and O is real and represents the imaginary
%           part of the transmittance, which is approximately considered as 
%           the targeted phase-shift image.
%
%       > 'absorbing': purely absorbing object. The linearized intensity
%       measurement formation model writes:
%
%           I = 1 + GZ * O
%
%           where GZ = 2 Re[HZ] and O is real and represents the opposite
%           of the opacity : O = T - 1.
%
%   According to the parameters, the kernel is based on an analytic
%   calculation either in the spatial or Fourier domain to respect the
%   Shannon's sampling theorem.
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
%


flag_phaseref = false;  % if not specified:
                        % default value for FLAG_PHASEREF
type_kernel = 'Hz';

if (nargin > 6)
    flag_phaseref = varargin{1};
    type_kernel = 'Hz';
    if (nargin > 8)
        % Check TYPE_HOLOGRAM and FLAG_LINEARIZE
        if (strcmp(varargin{2},'intensity') && varargin{3} == true)
            flag_phaseref = false;
            type_kernel = 'Gz_dephasing';
            % Check TYPE_OBJ
            if (nargin>9 && strcmp(varargin{4},'absorbing'))
                type_kernel = 'Gz_absorbing';
            end
        end
    end
end

if (strcmp(type_kernel,'Hz'))
    disp('Full Fresnel kernel Hz.');
elseif (strcmp(type_kernel,'Gz_dephasing'))
    disp('Adapted Gz kernel for reconstructing purely and weakly dephasing objects from intensity measurements.');
elseif (strcmp(type_kernel,'Gz_absorbing'))
    disp('Adapted Gz kernel for reconstructing purely absorbing objects from intensity measurements.');
end

% Wave number
k0 = 2*n_0*pi/lambda;

% Get z limit for Shannon sampling compatibility
z_shannon=max([npix_W, npix_H])*(pixel_size^2)/lambda;
fprintf('Shannon''s parameter = %.3f.\n', z_shannon);

flag_shannon = (z>=z_shannon);

if(flag_shannon)
    % Analytic impulse response
    disp('Impulse response.');
    
    x=get_fft_compatible_coordinates(npix_W, pixel_size);
    y=get_fft_compatible_coordinates(npix_H, pixel_size);
    [X,Y]=meshgrid(x,y);
    
    n0_on_i_pi_lambda_z = n_0/1i/lambda/z;
    hz=n0_on_i_pi_lambda_z*exp(n_0*1i*pi/lambda/z*(X.^2+Y.^2));
    if (flag_phaseref)
        hz=hz.*exp(1i*k0*z);
    end
    
    if (strcmp(type_kernel,'Gz_dephasing'))
        hz = -2*imag(hz);
    elseif (strcmp(type_kernel,'Gz_absorbing'))
        hz = 2*real(hz);
    end
    
    Hz=(pixel_size^2)*fft2(ifftshift(hz));
    
    if (~strcmp(type_kernel,'Hz'))
        Hz = real(Hz);
    end
    
    % Deprecated implementation (get the backpropagator)
    % H_z=conj(Hz);
else
    % Analytic transfert function (= Fourier transform of the impulse
    % response
    disp('Transfert function.');
    
    u=get_fft_compatible_coordinates(npix_W, 1.0/npix_W/pixel_size);
    v=get_fft_compatible_coordinates(npix_H, 1.0/npix_W/pixel_size);
    [U,V]=meshgrid(u,v);
    
    i_pi_lambda_z=-1i*pi*lambda*z/n_0;
    
    Hz=ifftshift(exp(i_pi_lambda_z*( U.^2+V.^2) ));
    
    if (flag_phaseref)
        Hz=Hz.*exp(1i*k0*z);
    end
    
    if (strcmp(type_kernel,'Gz_dephasing'))
        Hz = real(fft2(-2*imag(ifft2(Hz))));
    elseif (strcmp(type_kernel,'Gz_absorbing'))
        Hz = real(fft2(2*real(ifft2(Hz))));
    end
    
    % Deprecated implementation (get the backpropagator)
    % H_z=conj(Hz);
end

end
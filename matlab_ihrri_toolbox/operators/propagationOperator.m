function [xconv] = propagationOperator(x,H,varargin)
% [xconv] = propagationOperator(x,H,varargin)
%
%   This function performs a convolution of image X by the kernel H in the
%   Fourier space. 
%
%   In VARGARIN, 2 parameters can be given:
%   - FLAG_PAD: if true (default: true), a zero-padding will be performed 
%   provided that the size of H is doubled compared with the size of X.
%   - EXT_VAL: to set only if FLAG_TRANSP=false (default: 0). It specifies 
%   the value to set outside the initial field of view. 
%
%   Note: the kernel H is already "fftshifted".
%
% Created: 04/06/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516,
%           Université Jean Monnet,
%           F-42000 Saint-Étienne,
%           France
%           fabien.momey@univ-st-etienne.fr
%

flag_pad = true;
ext_val = 0;
if (nargin>2)
    flag_pad = varargin{1};
    if ((flag_pad) && (nargin > 3))
        ext_val = varargin{2};
    end
else
end

if (flag_pad)
    xconv = ifftshift(ifft2(fft2(ifftshift(fovExtensionOperator(x,2,0,ext_val))).*H));
    xconv = fovExtensionOperator(xconv,2,1);
else
    xconv = ifftshift(ifft2(fft2(ifftshift(x)).*H));
end

end

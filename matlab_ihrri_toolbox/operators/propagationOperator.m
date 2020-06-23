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

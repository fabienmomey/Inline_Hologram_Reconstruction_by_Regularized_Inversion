function [x_ext] = fovExtensionOperator(x,fov_extension_factor,varargin)
% [x_ext] = fovExtensionOperator(x,fov_extension_factor,varargin)
%
%   This function performs a zero-padded extension of the field of view
%   with a factor FOV_EXTENSION_FACTOR of image X.
%
%   In VARGARIN, 2 parameters can be given:
%   - FLAG_TRANSP: if true (default: false), the transpose operator is 
%   performed, i.e. the cropping operator. Then the image X has already a 
%   doubled size from the original field of view.
%   - EXT_VAL: to set only if FLAG_TRANSP=false (default: 0). It specifies 
%   the value to set outside the initial field of view. 
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

flag_transp = false;
ext_val = 0;
if (nargin>2)
    flag_transp = varargin{1};
    if ((~flag_transp) && (nargin > 3))
        ext_val = varargin{2};
    end
else
end

if (fov_extension_factor<1)
    
    error('An extension factor greater or equal to 1 must be given.');
    
elseif (fov_extension_factor>1)
    
    % Get image size
    [npix_H,npix_W] = size(x);
    
    if (~flag_transp)
        % Get extended size of the field of view
        % => integer part of fov_extension_factor*SIZE
        npix_H_ext_float = fov_extension_factor*npix_H;
        npix_W_ext_float = fov_extension_factor*npix_W;
        npix_H_ext = floor(npix_H_ext_float);
        npix_W_ext = floor(npix_W_ext_float);
        
        % Extended image (zero-padded)
        x_ext = ext_val*switch_to_gpu_array(ones(npix_H_ext,npix_W_ext));
        
        % Get starting pixel for filling the "useful" field of view
        start_H = floor((npix_H_ext-npix_H)/2.0);
        start_W = floor((npix_W_ext-npix_W)/2.0);
        
        if (mod(npix_H,2)==1 && mod(npix_H_ext,2)==0)      
            start_H = start_H + 1;         
        end
        
        if (mod(npix_W,2)==1 && mod(npix_W_ext,2)==0)      
            start_W = start_W + 1;         
        end
        
        x_ext(start_H+1:start_H+npix_H,...
                start_W+1:start_W+npix_W) = x;
    else
        npix_H_fov = ceil(npix_H/fov_extension_factor);
        npix_W_fov = ceil(npix_W/fov_extension_factor);
        
        % Get starting pixel for filling the "useful" field of view
        start_H = floor((npix_H-npix_H_fov)/2.0);
        start_W = floor((npix_W-npix_W_fov)/2.0);
        
        if (mod(npix_H_fov,2)==1 && mod(npix_H,2)==0)      
            start_H = start_H + 1;         
        end
        
        if (mod(npix_W_fov,2)==1 && mod(npix_W,2)==0)      
            start_W = start_W + 1;         
        end
        
        x_ext = x(start_H+1:start_H+npix_H_fov,...
                start_W+1:start_W+npix_W_fov);
    end
    
else
    
    x_ext = x ;
    
end

end
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
%

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
        x_ext = ext_val*ones(npix_H_ext,npix_W_ext);
        
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
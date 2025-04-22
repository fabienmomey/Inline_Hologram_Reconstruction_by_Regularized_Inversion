function [fx,gx,varargout] = critWLSlinearDenoising(x,y,M,M_,varargin)
% [fx,gx,varargout] = critWLSlinear(x,y,Gz,G_z,varargin)
%
%   This function computes the cost and gradient in X of a Weighted Least
%   Squares (WLS) criterion with a linear model:
%
%                       F(X)    = || C*M(X) - Y ||_W^2
%
%   The gradient Grad(X) writes:
%
%                                                 /            \         
%                       Grad(X)    = 2 * C * M^T  | C*M(X) - Y |
%                                                 \            /
%
%   X: current image guess depending on the type of object:
%           > 'dephasing': purely and weakly dephasing object.
%               => the unkown image X is real and represents the imaginary
%               part of the transmittance, which is approximately
%               considered as the targeted phase-shift image.
%               => the dimensions of X are just [width,height].
%           > 'absorbing': purely absorbing object.
%               => the unkown image X is real and represents the opposite
%               of the opacity : X = T - 1.
%               => the dimensions of X are just [width,height].
%
%   Y: data image (intensity measurements).
%
%   M: function handle to perform the forward operator = convolutive kernel 
%      (see propagationOperator functions).
%   M_: function handle to perform the backward operator = convolutive kernel 
%       (see propagationOperator functions).
%
%   In VARGARIN, 2 parameters can be given:
%   - Reg:  a function handle to a smooth regularizarion penalty
%   - C:    a constant scaling factor that accounts for the intensity of 
%           the incident wave |a_0|^2 as well as the detector gain and 
%           quantum efficiency [2]. If not set or <0, the scaling factor
%           can be computed "on-the-fly":
%
%           C = <M(X),Y> / <M(X),M(X)>
%
%           where <A,B> stands for the scalar product of vectors A and B.
%
%   - W:    diagonal elements of the inverse noise covariance matrix C^{-1}
%           => under hypothesis of uncorrelated noise [2].
%
%   The function returns:
%   - FX : the cost value (scalar)
%   - GX : the gradient image relative to X
%
%   In VARGAROUT, 2 additional parameters can be extracted:
%   - [fx,gx,c] = critWLS(x,y,M,M_,c) if c < 0 or c not set
%   - [fx,gx,residues] = critWLS(x,y,M,M_,c) if c > 0
%   - [fx,gx,c,residues] = critWLS(x,y,M,M_,c)
%
% Created: 04/15/2020 (mm/dd/yyyy)
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
[npix_W, npix_H] = size(x);

%% Is a smooth regularizer given ?
flag_reg = false ;
if (nargin>4 && ~isempty(varargin{1}))
    reg = varargin{1} ;
    flag_reg = true ;
else
    disp('No regularizer set.')
end

%% Criterion parameters
flag_c = true;
if (nargin>5)
    c = varargin{2};
    if (c > 0)
        flag_c = false;
    end
end

flag_w = false;
if (nargin>6)
    w = varargin{3};
    flag_w = true;
end

%% Forward propagation
uopt = M(x);
     
%% Calculate optimal scaling parameter c
if (flag_c && ~flag_w)
    c = sum(uopt(:).*y(:))/sum(uopt(:).*uopt(:)) ;
elseif (flag_c && flag_w)
    c = sum(uopt(:).*w(:).*y(:))/sum(uopt(:).*w(:).*uopt(:)) ;
end

%% Cost
cuopt_y = c*uopt-y;
if (flag_w)
    cuopt_y = w.*cuopt_y;
end
fx = sum(cuopt_y(:).^2);

%% Gradient
gx = 2*c*M_(cuopt_y);

%% Regul
if (flag_reg)
    [fxreg,gxreg] = reg(x) ;
    fx = fx + fxreg ;
    gx = gx + gxreg ;
end

if (nargout > 2)
    if (nargout > 3)
        varargout{1} = c;
        varargout{2} = cuopt_y;
    else
        if (flag_c)
            varargout{1} = c;
        else
            varargout{1} = cuopt_y;
        end
    end
end

end
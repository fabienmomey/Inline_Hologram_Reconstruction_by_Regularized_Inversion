function [xproj, fxproj] = softThresholdingOperator(x,mu,varargin)
% [xproj, fxproj] = softThresholdingOperator(x,mu,varargin)
%
%   This function performs a soft-thresholding operator on image X.
%
%   In VARGARIN, 3 parameters can be given:
%   - FLAG_CONST: -1, 0 or 1 (default: 0). If set to 1 (resp. -1), the
%   operator enforces positivity (resp. negativity) constraint on X.
%   - FLAG_CPLX: if true (default: false), the image X is complex-valued
%   the operator has to be applied on both the real and imaginary parts of
%   X, following the rule given by the next parameter FLAG_SEPARABLE. Note
%   the the parameter FLAG_CONST does not hold if FLAG_CPLX is true.
%   - FLAG_SEPARABLE: to set only if FLAG_CPLX=true (default: false). If
%   true, a soft-thresholding is applied separately to both the real and
%   imaginary part. This corresponds to the following regularizer:
%
%                                                          
%       Reg(X)  = mu * ( || Re[X] ||_1 + || Im[X] ||_1 )
%
%                            /                       \   
%               = mu * Sum_k | |Re[X_k]| + |Im[X_k]| |
%                            \                       / 
%
%   The operator's formula writes [1]:
%
%       Op(X_k) = sign(X_k) * max (0 , |X_k - mu|) 
%
%   
%   If false, the operator becomes the proximal operator of the euclidian
%   norm:
%
%       Reg(X)  = mu * || |X| ||_1 
%
%                            /                                 \   
%               = mu * Sum_k | sqrt(|Re[X_k]|^2 + |Im[X_k]|^2) |
%                            \                                 /  
%
%   which can be assimilated to the L1 norm of the modulus of X. The 
%   operator's formula writes [2]:
%
%                   /
%       Op(X_k) =   | (1 - mu/|X_k|) * X_k      if |X_k| > mu
%                   |
%                   | 0                         else
%                   \
%
% References
%
% - [1] A. Beck and M. Teboulle, “Fast gradient-based algorithms for con-
%       strained total variation image denoising and deblurring problems,”
%       IEEE Trans. Image Process. 18, 2419–2434, 2009.
% - [2] P.L.Combettes, J.-C. Pesquet, “A proximal decomposition method 
%       for solving convex variational inverse problems,” Inverse Probl., 
%       24(6):x+27, 2008.
%
% Created: 05/27/2020 (mm/dd/yyyy)
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

flag_const = 0.0;
flag_cplx = false;
if (nargin>2)
    flag_const = varargin{1};
    if (flag_const ~= 0.0)
        flag_const = flag_const/abs(flag_const);
    end
    if (nargin>3)
        flag_cplx = varargin{2};
        flag_separable = false;
        if (flag_cplx==true && nargin>4)
            flag_separable = varargin{3};
        end
    end
end

if (~flag_cplx) %% Case real
    if (flag_const==0.0)
        xproj = sign(x) .* max(0.0,abs(x)-mu);
    else
        xproj = flag_const * max(0.0,flag_const*(x-mu));
    end
    fxproj = mu*sum(abs(xproj(:))) ;
else            %% Case complex
    xproj = zeros(size(x));
    rex = x(:,:,1);
    imx = x(:,:,2);
    if (flag_separable)
        xproj(:,:,1) = sign(rex) .* max(0.0,abs(rex)-mu);
        xproj(:,:,2) = sign(imx) .* max(0.0,abs(imx)-mu);
        fxproj = mu*(sum(abs(xproj(:,:,1)),'all') + sum(abs(xproj(:,:,2)),'all')) ;
    else
        maxx = max(0.0,1.0-(mu./sqrt(rex.^2+imx.^2)));
        xproj(:,:,1) = maxx.*rex;
        xproj(:,:,2) = maxx.*imx;
        fxproj = mu*sum(sqrt(xproj(:,:,1).^2+xproj(:,:,2).^2),'all');
    end
end

end
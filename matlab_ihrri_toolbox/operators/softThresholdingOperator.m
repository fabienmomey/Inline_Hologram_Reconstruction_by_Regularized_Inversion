function [xproj] = softThresholdingOperator(x,mu,varargin)
% [xproj] = softThresholdingOperator(x,mu,varargin)
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
%

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
else            %% Case complex
    xproj = zeros(size(x));
    rex = x(:,:,1);
    imx = x(:,:,2);
    if (flag_separable)
        xproj(:,:,1) = sign(rex) .* max(0.0,abs(rex)-mu);
        xproj(:,:,2) = sign(imx) .* max(0.0,abs(imx)-mu);
    else
        maxx = max(0.0,1.0-(mu./sqrt(rex.^2+imx.^2)));
        xproj(:,:,1) = maxx.*rex;
        xproj(:,:,2) = maxx.*imx;
    end
end

end
function [fx,gx,varargout] = critWLSlinear(x,y,Gz,G_z,varargin)
% [fx,gx,varargout] = critWLSlinear(x,y,Gz,G_z,varargin)
%
%   This function computes the cost and gradient in X of a Weighted Least
%   Squares (WLS) criterion with a linear model:
%
%                       F(X)    = || C*M(X) - Y ||_W^2
%
%                               = || C*(1+G.X) - Y ||_W^2
%
%   The gradient G(X) writes:
%
%                                              /               \         
%                       G(X)    = 2 * C * G^T  | C*(1+G.X) - Y |
%                                              \               /
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
%   Y: data image (square root of intensity measurements).
%
%   Gz: function handle to perform the propagation operator  (see 
%       getFresnelPropagation and propagationOperator functions).
%   G_z: function handle to perform the backpropagation operator (see 
%       getFresnelPropagation and propagationOperator functions).
%
%   In VARGARIN, 2 parameters can be given:
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
%   In VARGAROUT, 2 additional parameters can be extracted:
%   - [fx,gx,c] = critWLS(x,y,Hz,H_z,c) if c < 0 or c not set
%   - [fx,gx,residues] = critWLS(x,y,Hz,H_z,c) if c > 0
%   - [fx,gx,c,residues] = critWLS(x,y,Hz,H_z,c)
%
% References
%
% - [1] J. Fienup, “Phase retrieval algorithms: a comparison,”
%                   Applied optics, vol. 21, no. 15, pp. 2758–2769, 1982.
% - [2] F. Momey, L. Denis, T. Olivier, C. Fournier, "From Fienup’s phase 
%                   retrieval techniques to regularized inversion for 
%                   in-line holography: tutorial," JOSA A, vol. 36, no. 12, 
%                   D62-D80, 2019. 
%
% Created: 06/12/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516,
%           Université Jean Monnet,
%           F-42000 Saint-Étienne,
%           France
%           fabien.momey@univ-st-etienne.fr
%

%% Extract size (in pixels) of the field of view
[npix_W, npix_H] = size(y);

%% Criterion parameters
flag_c = true;
if (nargin>4)
    c = varargin{1};
    if (c > 0)
        flag_c = false;
    end
end

flag_w = false;
if (nargin>5)
    w = varargin{2};
    flag_w = true;
end

%% Forward propagation
uopt = Gz(x)+1.0;
     
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
gx = 2*c*G_z(cuopt_y);

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
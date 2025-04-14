function [fx,gx,varargout] = critWLS(x,y,Hz,H_z,varargin)
% [fx,gx,varargout] = critWLS(x,y,Hz,H_z,varargin)
%
%   This function computes the cost and gradient in X of a Weighted Least
%   Squares (WLS) criterion with a nonlinear model of intensity
%   measurements:
%
%                  F(X)    = || C*M(X) - Y ||_W^2
%
%                          = || C*|1+H.X|^2 - Y ||_W^2
%
%   The gradient G(X) writes:
%                                        /                               \
%                                        |           /                 \ |
%                  G(X)    = 4 * C * H^T | (1+H.X) . | C*|1+H.X|^2 - Y | |
%                                        |           \                 / |
%                                        \                               /
%
%   X: current image guess => a 2-component image ([width,height,2]), each
%   one corresponding respectively to the real and imaginary part of the
%   complex deviation from the unit transmittance plane)
%
%   Y: data image (intensity measurements).
%
%   Hz: function handle to perform the propagation operator  (see 
%       getFresnelPropagation and propagationOperator functions).
%   H_z: function handle to perform the backpropagation operator (see 
%       getFresnelPropagation and propagationOperator functions).
%
%   In VARGARIN, 3 parameters can be given:
%   - C:    a constant scaling factor that accounts for the intensity of 
%           the incident wave |a_0|^2 as well as the detector gain and 
%           quantum efficiency [2]. If not set or <0, the scaling factor
%           can be computed "on-the-fly" (and returned in VARARGOUT):
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
% Created: 06/10/2020 (mm/dd/yyyy)
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
[npix_W, npix_H, ncomp] = size(x);

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
uopt = Hz(x(:,:,1)+1i*x(:,:,2))+1.0;
Iopt = abs(uopt).^2;
     
%% Calculate optimal scaling parameter c
if (flag_c && ~flag_w)
    c = sum(Iopt(:).*y(:))/sum(Iopt(:).*Iopt(:)) ;
elseif (flag_c && flag_w)
    c = sum(Iopt(:).*w(:).*y(:))/sum(Iopt(:).*w(:).*Iopt(:)) ;
end

%% Cost
cIopt_y = c*Iopt-y;
if (flag_w)
    cIopt_y = w.*cIopt_y;
end
fx = sum(cIopt_y(:).^2);

%% Gradient
gxcplx = 4*c*H_z((uopt.*(cIopt_y)));

gx = zeros(npix_W,npix_H,2);
%gx(:,:,1) = real(4*c*H_z(real(uopt).*(cIopt_y)));
%gx(:,:,2) = imag(4*c*H_z(imag(uopt).*(cIopt_y)));

gx(:,:,1) = real(gxcplx);
gx(:,:,2) = imag(gxcplx);

if (nargout > 2)
    if (nargout > 3)
        varargout{1} = c;
        varargout{2} = cIopt_y;
    else
        if (flag_c)
            varargout{1} = c;
        else
            varargout{1} = cIopt_y;
        end
    end
end

end
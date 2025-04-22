function [fx,gx] = critEdgePreservingSmoothing(x,Gx,Gy,GxT,GyT,mu,epsilon,varargin)
% [fx,gx] = critEdgePreservingSmoothing(x,G,varargin)
%
%   This function computes ...
%
%   X: current image guess => a 2-component image ([width,height,2]), each
%   one corresponding respectively to the real and imaginary part of the
%   complex deviation from the unit transmittance plane)
%
%   Gx: function handle to perform the gradient operator in x direction.
%   Gy: function handle to perform the gradient operator in y direction.
%   GxT: function handle to perform the adjoint gradient operator in
%        x direction.
%   GyT: function handle to perform the adjoint gradient operator in
%        y direction.
%   (see getGradientOperator and propagationOperator functions)
%
%   In VARGARIN, 2 parameters can be given:
%   - FLAG_CPLX: if true (default: false), the image X is complex-valued
%   the operator has to be applied on both the real and imaginary parts of
%   X, following the rule given by the next parameter FLAG_SEPARABLE. Note
%   the the parameter FLAG_CONST does not hold if FLAG_CPLX is true.
%   - FLAG_SEPARABLE: to set only if FLAG_CPLX=true (default: false). If
%   true, the edge preserving regularizer is applied separately
%   to both the real and imaginary part.
%   This corresponds to the following regularizer:
%
%
%
%   The function returns:
%   - FX : the cost value (scalar)
%   - GX : the gradient image relative to X
%
% Created: 11/04/2025 (mm/dd/yyyy)
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

flag_cplx = false;
if (nargin>7)
    flag_cplx = varargin{1};
    flag_separable = false;
    %% Weigth hyperparameter (can be different for the parts)
    if (~isscalar(mu) && (all(size(mu)==[1,2]) || all(size(mu)==[2,1])))
        mu_real = mu(1) ;
        mu_imag = mu(2) ;
        flag_separable = true ;
    else % else : keep using mu
        mu = mu(1) ;
    end
    %% Threshold hyperparameter (can be different for the parts)
    if (~isscalar(epsilon) && flag_separable && (all(size(epsilon)==[1,2]) || all(size(epsilon)==[2,1])))
        epsilon_real = epsilon(1) ;
        epsilon_imag = epsilon(2) ;
    else % else : keep using epsilon
        epsilon = epsilon(1) ;
    end
end

if (~flag_cplx) %% Case real
    %% Calculate the gradient components
    grad_x = Gx(x) ;
    grad_y = Gy(x) ;
    %% Compute the cost
    cost_image = sqrt(grad_x.^2+grad_y.^2+epsilon^2) ;
    fx = mu * sum(cost_image(:)) ;
    %% Compute the gradient
    gx = mu*(GxT(grad_x)+GyT(grad_y))./cost_image ;
else  %% Case complex
    %% Calculate the gradient components
    grad_x_real = Gx(x(:,:,1)) ;
    grad_y_real = Gy(x(:,:,1)) ;
    grad_x_imag = Gx(x(:,:,2)) ;
    grad_y_imag = Gy(x(:,:,2)) ;
    %% Compute the cost and gradient
    if (flag_separable)
        %% Compute the cost
        cost_image_real = sqrt(grad_x_real.^2+grad_y_real.^2+epsilon_real^2) ;
        cost_image_imag = sqrt(grad_x_imag.^2+grad_y_imag.^2+epsilon_imag^2);
        fx = mu_real*sum(cost_image_real(:)) + mu_imag*sum(cost_image_imag(:)) ;
        %% Compute the gradient
        gx = mu_real*(GxT(grad_x_real)+GyT(grad_y_real))./cost_image_real + ...
            mu_imag*(GxT(grad_x_imag)+GyT(grad_y_imag))./cost_image_imag ;
    else
        %% Compute the cost
        cost_image = sqrt(grad_x_real.^2 + grad_y_real.^2 + ...
            grad_x_imag.^2 + grad_y_imag.^2 + epsilon^2) ;
        fx = mu * sum(cost_image(:)) ;
        %% Compute the gradient
        gx = mu*(GxT(grad_x_real)+GyT(grad_y_real)+GxT(grad_x_imag)+GyT(grad_y_imag))./cost_image ;
    end
end

end
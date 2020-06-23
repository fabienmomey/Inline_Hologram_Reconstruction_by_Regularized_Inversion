% Test script of function:
%
% [x_ext] = fovExtensionOperator(x,fov_extension_factor,flag_transp)
%
% This script uses as a meshgrid of 2D spatial coordinates with the (0,0) 
% coordinate located exactly on a pixel position.
%
% The test helps verifying that an application of a direct operator
% (padding operator) followed by its adjoint (cropping operator) is
% equivalent to perform the identity operator.
%
% Moreover, the function has been implemented so that performing FFTSHIFT
% (or IFFTSHIFT in the case of odd dimensions) moves the (0,0) coordinate 
% at the top-left corner of the meshgrid matrix to ensure coherent
% application of FFT algorithms.
%
% Created: 05/12/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516,
%           Université Jean Monnet,
%           F-42000 Saint-Étienne,
%           France
%           fabien.momey@univ-st-etienne.fr
%
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

clear all;
close all;
clc;

addpath(genpath('../.'));

%% Choose test parameters
extension_factor = 2.3;
npix_W = 5;
npix_H = 4;

x = (0:npix_W-1)-floor(npix_W/2.0);
y = (0:npix_H-1)-floor(npix_H/2.0);

[X,Y] = meshgrid(x,y);

%% Padding operator
Xext = fovExtensionOperator(X,extension_factor,0,-1);
Yext = fovExtensionOperator(Y,extension_factor,0,-1);

%% Cropping operator
Xfov = fovExtensionOperator(Xext,extension_factor,1);
Yfov = fovExtensionOperator(Yext,extension_factor,1);
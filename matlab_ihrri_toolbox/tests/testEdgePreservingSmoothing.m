% Test script of function:
%
% [fx,gx] = critEdgePreservingSmoothing(x,Gx,Gy,GxT,GyT,mu,epsilon,varargin)
%
% on a simple denoising task
%
% Created: 04/15/2025 (mm/dd/yyyy)
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

%% Create the dataset
radius_o = 64.0 ;
npix_W = 512 ;
npix_H = 512 ;
pixel_size = 1.0 ;
x=get_fft_compatible_coordinates(npix_W, pixel_size);
y=get_fft_compatible_coordinates(npix_H, pixel_size);
[X,Y]=meshgrid(x,y);

m = fspecial('gaussian',[npix_H*2,npix_W*2],2) ;
M = fft2(ifftshift(m)) ;

data_unnoisy = double(X.^2+Y.^2<radius_o^2) ;
ihrri_show(data_unnoisy,'Data') ;

data_unnoisy = propagationOperator(data_unnoisy,M,true) ;

sigma_noise = 0.1 ;
data_noisy = data_unnoisy + sigma_noise*randn(npix_H,npix_W) ;

ihrri_show(data_noisy,'Data') ;

%% Forward and backward models
Forward = @(o) (propagationOperator(o,M,true));
Backward = @(o) (propagationOperator(o,conj(M),true));

%% Regularization
muEdgePres = 0.001 ;
epsilonEdgePres = 1.0e-4 ;
Grad = getGradientOperator(2*npix_W, 2*npix_H) ;
% Gradient in X direction
Gx = @(o) (propagationOperator(o,Grad,true)) ;
% Gradient in Y direction
Gy = @(o) (propagationOperator(o,conj(Grad'),true)) ;
% Adjoint gradient in X direction
GxT = @(o) (propagationOperator(o,conj(Grad),true)) ;
% Adjoint gradient in Y direction
GyT = @(o) (propagationOperator(o,Grad',true)) ;
Reg = @(x) (critEdgePreservingSmoothing(x,Gx,Gy,GxT,GyT,muEdgePres,epsilonEdgePres)) ;

% Vérif. correct "adjointness" of gradient operators
% x=randn(npix_H, npix_W) ;
% y=randn(npix_H, npix_W) ;
% sum(sum(y.*Gx(x)))
% sum(sum(x.*GxT(y)))
% sum(sum(y.*Gy(x)))
% sum(sum(x.*GyT(y)))

%% Global criterion
Crit = @(x) (critWLSlinearDenoising(x,data_noisy,Forward,Backward,Reg,1)) ;

%% Reconstructor
RECoptions = struct('type_obj','absorbing',...
            'flag_linearize',true,...
            'rconst',[0,Inf],...
            'maxiter', 50,...
            'flag_cost',true,...
            'flag_evolx',false,...
            'verbose',true);

x0 = zeros(npix_H,npix_W) ;
[xopt,fxopt,Gxopt,c,evolcost] = algoRI(x0,Crit,RECoptions);
ihrri_show(xopt,'Reconstruction');
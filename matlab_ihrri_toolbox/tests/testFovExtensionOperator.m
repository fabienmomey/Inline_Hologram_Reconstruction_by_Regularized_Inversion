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
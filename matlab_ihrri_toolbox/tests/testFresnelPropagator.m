% Test script of functions:
%
% [Hz] = getFresnelPropagator(npix_W, npix_H, pixel_size, lambda, n_0, z, varargin)
% [xconv] = propagationOperator(x,H,varargin)
%
% This script loads a ground truth image (deviation from the 1-transmittance plane) 
% and perform the propagation operator with and without 0-padding. Then it
% applies the backpropagation of the total diffracted wave UTOT.
%
% Created: 05/21/2020 (mm/dd/yyyy)
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

datadir = '../data/Bead_simulations/';

% Load image
run([datadir,'parameters.m']);
o = fitsread([datadir,'ground_truth_GaussBeads_Sept2019_pure_phase_FR_mag5.67e+01_z1.25e-05.fits']);
o = o(:,:,1) + 1i*o(:,:,2);
%o = fovExtensionOperator(o,2,true);
osize = size(o);

% Parameters
z_s = 7.2822e-6 ;
mag = 56.7; 
lambda = 532e-9;
n_0 = 1.52; 
pixel_size = 2.2e-6 / mag; 
npix_H = osize(1); 
npix_W = osize(2); 

x = get_fft_compatible_coordinates(npix_W, pixel_size);
y = get_fft_compatible_coordinates(npix_H, pixel_size);
[X,Y] = meshgrid(x,y);

ihrri_show(imag(o), 'Phase Beads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Propagation test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Without zero-padding

[Hz]=getFresnelPropagator(npix_W, npix_H, pixel_size, lambda, n_0, z_s);

ihrri_show(fftshift(real(ifft2(Hz))), 'No padding: Real part impulse response');
ihrri_show(fftshift(imag(ifft2(Hz))), 'No padding: Imag part impulse response');

utot_nopad = 1.0+propagationOperator(o,Hz,false);

ihrri_show(abs(utot_nopad).^2, 'Intensity');

%% With zero-padding (the best way !)

[Hz_pad]=getFresnelPropagator(2*npix_W, 2*npix_H, pixel_size, lambda, n_0, z_s);

ihrri_show(fftshift(real(ifft2(Hz_pad))), '0-padding: Real part impulse response');
ihrri_show(fftshift(imag(ifft2(Hz_pad))), '0-padding: Imag part impulse response');

utot_pad = 1.0+propagationOperator(o,Hz_pad,true);

ihrri_show(abs(utot_pad).^2, '0-padding: Intensity');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Backpropagation test %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% With zero-padding (the best way !)

orec = propagationOperator(utot_pad-1.0,conj(Hz_pad),true);

ihrri_show(imag(orec), 'Backpropagation');



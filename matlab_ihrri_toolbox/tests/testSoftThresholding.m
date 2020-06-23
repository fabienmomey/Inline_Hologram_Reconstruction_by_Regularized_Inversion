% Test script of function:
%
% [xproj] = softThresholdingOperator(x,mu,varargin)
%
% Created: 05/27/2020 (mm/dd/yyyy)
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

%% Choose context parameters
flag_const = 0;
flag_cplx = true;
flag_separable = false;

%% Generate a random starting point
if (flag_cplx)
    x0 = randn(1,1) + 1i * randn(1,1);
    x0 = x0/abs(x0);
else
    x0 = 0.5+randn(1,1);
end

%% Generate a vector of hyperparameters mu
mu = 0:0.1:2;
nmu = length(mu);

%% Mapping criterion space
N = 512;
x = 4*(0:N-1)/(N-1)-2;

subplotsize = ceil(sqrt(nmu));
figure(1); hold on;
if (flag_cplx)
    [X,Y] = meshgrid(x,x);
    for k=1:nmu
        Crit = 0.5*((X-real(x0)).^2 + (Y-imag(x0)).^2);
        if (~flag_separable)
            Crit = Crit + mu(k)*sqrt(X.^2+Y.^2);
        else
            Crit = Crit + mu(k)*(abs(X) + abs(Y));
        end
        subplot(subplotsize,subplotsize,k); contourf(X,Y,Crit,20); axis xy; hold on;
        plot(real(x0),imag(x0),'og');
        %% Apply soft-thresholding
        xproj=softThresholdingOperator(x0,mu(k),flag_const,flag_cplx,flag_separable);
        plot(real(xproj),imag(xproj),'*r');
        clear Crit xproj;
    end
else
    for k=1:nmu
        Crit = 0.5*(x-x0).^2 + mu(k)*abs(x);
        maxCrit = max(Crit(:));
        subplot(subplotsize,subplotsize,k); plot(x,Crit,'-b'); hold on;
        plot([x0,x0],[0,maxCrit],'--g');
        %% Apply soft-thresholding
        xproj=softThresholdingOperator(x0,mu(k),flag_const);
        plot([xproj,xproj],[0,maxCrit],'--r');
        clear Crit xproj;
    end
end



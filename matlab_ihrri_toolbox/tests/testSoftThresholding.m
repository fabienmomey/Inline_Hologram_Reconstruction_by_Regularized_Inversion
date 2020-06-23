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



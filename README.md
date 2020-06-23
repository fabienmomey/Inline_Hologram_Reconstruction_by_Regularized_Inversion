# Inline Hologram Reconstruction by Regularized Inversion (IHRRI) toolbox
A matlab code for "inverse problems" based image reconstruction of digital in-line holograms. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This toolbox implements Inverse Problems based algorithms dedicated to
% image reconstruction in digital in-line holographic microscopy (DHM).
% Theoretical aspects on DHM and inverse approaches are developed in a
% tutorial published in JOSA A [1], and this code constitutes a 
% demonstrator of the algorithms presented in this publication.
%
% Here is the main reconstruction script that can be executed as is, and
% can perform a reconstruction from in-line hologram data whose parameters
% (data and results saving paths, calibration, algorithm settings) have to
% be set in the file "parameters.m" (refer to it for more details).
%
% Summary
% -------
%
% The code is able to perform 2 kind of "Inverse problems" algorithms
% aiming at reconstructing an image X from an intensity in-line hologram 
% image Y. In this code, X is 2-component image ([width,height,2]), each
% one corresponding respectively to the real and imaginary part of the
% complex deviation from the unit transmittance plane :
%
%                   T = 1 + X(:,:,1) + i X(:,:,2)
%
% The 2 implemented reconstruction algorithm are:
%   - The Fienup's Error-Reduction algorithm [2], which stands for a
%   gradient descent algorithm aiming at solving the following
%   problem:
%   
%   X   = ARG MIN   || C*M(X) - Y ||_2^2        s.t.    X in Omega
%            x
%
%   where 
%       - M(X) = |1 + H.X| is the square root intensity hologram formation 
%       model of X using a convolutive propagation operator H (the Fresnel 
%       kernel is available in this toolbox).
%       - Omega stands for the validity domain of X that is enforced as a
%       "projection on constraints" operator in the algorithm. The domain
%       Omega takes the form of bound constraints applied on the real and
%       imaginary parts of X:
%           * RCONST = [XRMIN,XRMAX]
%           * ICONST = [XIMIN,XIMAX]
%
%       - C is a scaling factor that accounts for the intensity of the 
%       incident wave |a_0|^2 as well as the detector gain and quantum 
%       efficiency [1].
%
%   - The FISTA algorithm [3], which is a proximal gradient descent
%   algorithm aiming at solving the following sparsity problem:
%   
%   X   = ARG MIN   || C*M(X) - Y ||_W^2 + mu * || X ||_1
%            x
%
%                                          s.t.    X in Omega
%
%   where 
%       - M(X) is still the intensity hologram formation model which can 
%       take 2 forms:
%           * Considering purely and weakly dephasing or purely absorbing 
%           objects, we can use a linearized intensity hologram formation
%           model:
%
%           M(X) = 1 + G.X ~ |1 + H.X|^2 
%
%           where X becomes purely real image (size [width,height]) and G
%           is a purely real kernel (see the code's documentation for
%           details) which depends on the convolutive propagation operator 
%           H.
%
%           * Considering a unknown object:
%
%           M(X) = |1 + H.X|^2 
%           
%           that models the "full" intensity hologram image from the
%           "complex" image X.
%
%       - Omega stands for the validity domain of X that is enforced as a
%       "projection on constraints" operator in the algorithm.
%       - C is a scaling factor that accounts for the intensity of the 
%       incident wave |a_0|^2 as well as the detector gain and quantum 
%       efficiency [1].
%       - W is the inverse noise covariance matrix C^{-1}.
%
%       - || X ||_1 is a regularizer enforcing a sparsity constraint,
%       weighted by a parameter MU, and applied as a soft-thresholding
%       operator in the algorithm.
%
% Note that the code is able to perform fied-of-view extension in a
% rigorous way by setting an extension factor in the settings (see
% "parameters.m").
%
% In the future, an edge-preserving regularizer will be implemented which
% can be added to the reconstruction criterion (see [1] for theoretical
% details).
%
% References
%
% - [1] F. Momey, L. Denis, T. Olivier, C. Fournier, "From Fienup’s phase 
%                   retrieval techniques to regularized inversion for 
%                   in-line holography: tutorial," JOSA A, vol. 36, no. 12, 
%                   D62-D80, 2019. 
% - [2] J. Fienup, “Phase retrieval algorithms: a comparison,”
%                   Applied optics, vol. 21, no. 15, pp. 2758–2769, 1982.
% - [3]  A. Beck and M. Teboulle, “Fast gradient-based algorithms for con-
%       strained total variation image denoising and deblurring problems,”
%       IEEE Trans. Image Process. 18, 2419–2434, 2009.
%
% Created: 04/06/2020 (mm/dd/yyyy)
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

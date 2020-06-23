function [x] = get_fft_compatible_coordinates(npix, pixel_size)
% [x] = get_fft_compatible_coordinates(npix, pixel_size)
%
%   This function calculates a vector X of pixel coordinates centered at 0 
%   given the number of pixels NPIX and the pixel size PIXEL_SIZE. 
%   The coordinates are arranged so that the 0 falls at an exact pixel 
%   position, and also that a FFTSHIFT or IFFTSHIFT leads this 0 to be the 
%   first element of the vector X.
%
% Created: 04/06/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516,
%           Université Jean Monnet,
%           F-42000 Saint-Étienne,
%           France
%           fabien.momey@univ-st-etienne.fr
%

x=((0:npix-1)-ceil(0.5*(npix-1)))*pixel_size;

end
function [] = ihrri_show(im, varargin)
% [] =  ihrri_show(im, varargin)
%
% "Hand-made" function for showing an image with the screen resolution.
%
% VARARGIN arguments:
%
%   - varargin{1} : figure title.
%
% See also figure, imagesc.
%
% Created: 05/21/2020 (mm/dd/yyyy)
% Author:   Fabien Momey
%           Laboratoire Hubert Curien UMR CNRS 5516,
%           Université Jean Monnet,
%           F-42000 Saint-Étienne,
%           France
%           fabien.momey@univ-st-etienne.fr
%

figtitle = ['IHRRI Figure'];
if (nargin <= 2)
    if (nargin==2)
        if (ischar(varargin{1}))
            figtitle = varargin{1};
        else
            warning('Wrong type for FIGURE TITLE: must be CHAR. Keep default name.');
        end
    end
else
    error('Wrong number of VARARGIN arguments: only the FIGURE TITLE (CHAR) can be specified');
end

% get screen size
scrsz = get(groot,'ScreenSize');
% get image size
imsize = size(im);

h = figure('Name',figtitle,'Position',[scrsz(3)/4 scrsz(4)/2 imsize(2) imsize(1)]);
imagesc(im); colormap(gray);
axis xy;
ax = gca;
ax.BoxStyle = 'full';
ax.Box = 'on';
ax.XTick = [];
ax.YTick = [];
set(ax,'Position',get(ax,'OuterPosition'));
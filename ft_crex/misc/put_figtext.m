function gt = put_figtext(sstr, loc, sz, txcol, bkcol)
% Put text on a current figure axis according to loc indication (see inputs)
%
% ------------
% Inputs:
%
% - sstr [REQUIRED]: text to plot (string) - if numerical value, num2str is applied on sstr 
%
% - loc: 
%       -- position indication as 2 cardinal points from the axis center
%           'nw' or 'wn' : NorthWest - top-left corner [default]
%           'ne' or 'en' : NorthEast - top-right corner
%           'sw' or 'ws' : SouthWest - bottom-left corner
%           'se' or 'es' : SouthEast - bottom-right corner
%           loc can be indicated too with the combinaison of 2 letters from 
%           't' (top), 'b' (bottom), 'l' (left) and 'r'(right)
%       -- or a vector of coordinates [x, y] depending on data 
%
% - sz : text size [default: 12]
% 
% - txcol : text color [default: [0 0 0]] (can be standard letter ex: 'b' for blue)
%
% - bkcol : background color [default: [1 1 1]]
%
%-CREx-180612

if nargin < 1 || isempty(sstr)
    return;
end

if isnumeric(sstr)
    sstr = num2str(sstr);
end

% Check for text color property
% Set default
if nargin < 5 || isempty(bkcol) || (~ischar(bkcol) && length(bkcol)<3)
    bkcol = [1 1 1];
end
% Set orange color for 'o' character
if ischar(bkcol) && strcmp(bkcol, 'o')
    bkcol = [0.8 0.3 0.1];
end

% Check for background color
% Set default
if nargin < 4 || isempty(txcol) || (~ischar(txcol) && length(txcol)<3)
    txcol = [0 0 0];
end
% Set orange
if ischar(txcol) && strcmp(txcol, 'o')
    txcol = [0.8 0.3 0.1];
end

% Check for text size property
% Set default
if nargin < 3 || isempty(sz)
    sz = 12;
end

% Check for text location on gca
% All posible locations
alloc = {'nw' 'ne' 'sw' 'se' 'tl' 'tr', 'bl', 'br'};
alloc = [alloc , cellfun(@(x) x([2 1]), alloc, 'uniformoutput', 0)];
if nargin < 2 || (ischar(loc) && sum(strcmpi(loc, alloc))==0)
    loc = 'nw';    
end

if ~ischar(loc)
    cx = loc(1);
    cy = loc(2);
else
    [cx, cy, uini] = det_pos(loc);
end

gt = text(cx, cy, sstr,...
    'HorizontalAlignment','left',...
    'verticalalignment','top',...
    'color',txcol,...
    'BackgroundColor',bkcol,...
    'fontsize',sz); 

if ischar(loc)
    % Add a small amount for text in bottom
    if containss(loc, 's') || containss(loc, 'b')
        ext = get(gt, 'extent');
        yl = ylim;
        if ext(2) < yl(1)
            cy = cy + diff(ext([2 4]));
        end
        set(gt, 'position', [cx cy + ext(end)])
    end
    if containss(loc, 'e') || containss(loc, 'r')
        set(gt, 'horizontalalignment','right')
    end

    set(gca, 'units', uini);
end

function [cx, cy, uini] = det_pos(loc)
loc = lower(loc);

xl = xlim;
Dx = diff(xl);
    
yl = ylim;
Dy = diff(yl);

dtk = get(gca,'ticklength');

uini = get(gca, 'units');


set(gca,'units','centimeter')
pos = get(gca,'position');

L = pos(3); 
H = pos(4);

% Set minimum margin around text based on tick length
if L >= H
    htk = Dx*dtk(1);
    vtk = (dtk(1)*L*Dy)/H;
else
    vtk = Dy*dtk(1);
    htk = (dtk(1)*H*Dx)/L;
end
dtx = 1.5 * htk;
dty = 1.5 * vtk;

isylog = strcmp(get(gca, 'YScale'), 'log');
isxlog = strcmp(get(gca, 'XScale'), 'log');
% Set text position from the top (y coordinate)
if containss(loc, 'n') || containss(loc, 't')
    if isylog
        yt = get(gca, 'YTick');
        cy = yl(2) - diff(yt(end-1:end))/10;
    else
        cy = yl(2) - dty;
    end
else
    if isylog
        yt = get(gca, 'YTick');
        cy = yl(1) + diff(yt(1:2))/10;
    else
        cy = yl(1) + dty;
    end
end
% Set text position from the left (x coordinate)
if containss(loc, 'w') || containss(loc, 'l')
    if isxlog
        xt = get(gca, 'XTick');
        cx = xl(1) + diff(xt(1:2))/10;
    else
        cx = xl(1) + dtx;
    end
else
    if isxlog
        xt = get(gca, 'XTick');
        cx = xl(2) - diff(xt(end-1:end))/10;
    else    
        cx = xl(2) - dtx;
    end
end

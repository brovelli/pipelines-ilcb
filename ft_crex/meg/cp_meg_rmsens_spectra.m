function rms = cp_meg_rmsens_spectra(pfig, rmsini)
% Get the bad sensors from interactive selection on the spectra figure
%-CREx180720

if nargin < 2
    rmsini = [];
else
    if ~isempty(rmsini)
        if ~iscell(rmsini)
            rmsini = {rmsini};
        end
        sz = size(rmsini);
        if sz(2) > 1
            rmsini = rmsini';
        end
    end
end
if ~exist(pfig, 'file')
    warning('%s to select bad channel doesn''t exist', pfig)
    rms = rmsini;
    return;
end

open(pfig);
init_col;
hlist = add_controls(rmsini);
uiwait();

rms = get_add(hlist);

delete(gcf)

% Add callbacks
function hlist = add_controls(rmsini)
hlist = find_annot('lab_bad');

% Add other channels that were declare as bad [outside the cleanup figure step]
hlist = add_rmsens(hlist, rmsini);

% Add validate button callback function
hval = findobj(gca, 'tag', 'but_val');
% Callback function when validate button is clicked
set(hval, 'ButtonDownFcn', {@validate});
% set(hval, 'ButtonDownFcn', {@validate});

% Add reset button callback function
hbad = findobj(gca, 'tag', 'but_bad');
% Callback function when validate button is clicked
set(hbad, 'ButtonDownFcn', {@add_bad, hlist});

% Add restore callback function
hgood = findobj(gca, 'tag', 'but_good');
set(hgood, 'ButtonDownFcn', {@rm_bad, hlist});

set(gcf, 'WindowButtonDownFcn', {@dispnamesp, hlist});

set(gcf, 'CloseRequestFcn', []);

function keep_bad(hlist)
rma = get_add(hlist);
Na = length(rma);
for i = 1 : Na
    cline = findobj(gca, 'displayname', rma{i});
    if ~isempty(cline)
        cline.Color = [0.75 0.73 0.73];
        cline.LineWidth = 1.1;
    end
end

function hlist = add_rmsens(hlist, rms)
if ~isempty(rms)
    hlist.String = unique([get_add(hlist); rms]);   
end
keep_bad(hlist);

function rmall = get_add(hlist)
badlab = hlist.String;
if isempty(badlab)
    rmall = [];
else
    rmall = badlab;
end

function add_bad(~, ~, hlist)  
ff = get(gcf, 'WindowButtonDownFcn');
cline = ff{2};
clab = get(cline, 'displayname');
cline = findobj('DisplayName', clab);
rmsa = get_add(hlist);
if ~isempty(clab)
    rmsa = unique([rmsa; {clab}]);
    hlist.String = rmsa;
    cline(1).Color = [0.75 0.73 0.73];
    cline(1).LineWidth = 1.1;
end
set(gcf, 'WindowButtonDownFcn', {@dispnamesp, hlist});

function rm_bad(~, ~, hlist)  
ff = get(gcf, 'WindowButtonDownFcn');
cline = ff{2};
clab = get(cline, 'displayname');
cline = findobj('DisplayName', clab);
rmsa = get_add(hlist);
if ~isempty(clab)
    rmsa = rmsa(~strcmp(rmsa, clab));
    hlist.String = rmsa;
    set(cline, cline.UserData{1}, cline.UserData{2});
end
set(gcf, 'WindowButtonDownFcn', {@dispnamesp, hlist});

% Validate
function validate(~, ~) 
uiresume;

function han = find_annot(tagnam)
han = findall(0, 'tag', tagnam);
Nh = length(han);
if Nh > 1
    for i = 1 : Nh
        if han(i).Parent.Parent == gcf
            han = han(i);
            break;
        end
    end
end

% Display channel name on figure title
function dispnamesp(~, ~, hlist)

keep_bad(hlist);
av = findobj(gca, 'type', 'line', 'color', [0 1 0.6]);
if ~isempty(av)
    set(av, av.UserData{1}, av.UserData{2})
end
% If object is of type line (marker), its size is slightly increased when
% it is selected and color changed to green-blue
gobj = gco;
if gobj==gca || gobj==gcf || (isprop(gobj, 'type') && strcmp(gobj.Type, 'text'))
    title('');
end
if isprop(gobj, 'type') && isprop(gobj, 'displayname') && ~isempty(get(gobj, 'displayname'))

    dnam = get(gobj, 'displayname');

    % Keep initial properties and set new properties to highlight object
    iniprop = {{'linewidth', 'color'}, {get(gobj,'linewidth'), get(gco, 'color')}};
    curprop = {{'linewidth', 'color'}, {get(gobj,'linewidth')+0.2, [0 1 0.6]}};
    
    title(dnam,'fontsize',24,'interpreter','none');
    set_prop(gobj, curprop);

    set(gcf, 'WindowButtonDownFcn', {@after_click, gobj, iniprop});
    set(gcf, 'WindowButtonUpFcn', {@dispnamesp, hlist});
end

function after_click(~, ~, gobj, prop)
set_prop(gobj, prop)
   
function set_prop(gobj, prop)
Np = length(prop(:,1));
for i = 1 : Np
    set(gobj, prop{i, 1}, prop{i, 2});
end

function init_col
adisp = findobj(gca, 'type', 'line');
adisp = adisp(isprop(adisp, 'displayname'));

Na = length(adisp);
for i = 1 : Na
    set(adisp(i), 'UserData', {{'Color', 'LineWidth'}, {adisp(i).Color, adisp(i).LineWidth}});
end
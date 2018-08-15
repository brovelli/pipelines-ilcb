function dispname(~, ~)
% Function to display name of objects in a figure when the object is selected with
% the mouse (objects with 'displayname' property set to the name of the
% object)
% Name appears as figure title.
% In the main script where figure is created the command:
% set(gcf, 'WindowButtonDownFcn', @dispname);
% is added to allow the function dispname to be excecute each time an
% object is selectionned.
% 
%-CREx170118-- ft_CREx toolbox

% first = isfirst;

% Background color to define color of the title
if all(get(gcf,'color'))
    titcol = [0 0 0];
else
    titcol = [1 1 1];
end

% If object is of type line (marker), its size is slightly increased when
% it is selected and color changed to green-blue
gobj = gco;
if gobj==gca || gobj==gcf || (isprop(gobj, 'type') && strcmp(gobj.Type, 'text'))
    title('');
end
if isprop(gobj, 'type') && isprop(gobj, 'displayname') && ~isempty(get(gobj, 'displayname'))
    typ = get(gobj, 'type');
    dnam = get(gobj, 'displayname');

    % Keep initial properties and set new properties to highlight object
    if strcmp(typ, 'line')
        iniprop = {'markersize', get(gobj,'markersize')
                'color', get(gco, 'color')};
        curprop = {'markersize', get(gobj,'markersize')+2
                'color', [0 1 0.6]};

    % If object is a patch, facealpha is set to 1 and color to green-blue
    elseif strcmp(typ, 'patch')
        iniprop = {'facealpha', get(gobj,'facealpha')
                'facecolor', get(gobj, 'facecolor')};
        curprop = {'facealpha', 0.3
                'facecolor', [0 1 0.6]};    
    end

    
    title(dnam,'fontsize',24, 'color', titcol,'interpreter','none');

    
    set_prop(gobj, curprop);

    set(gcf, 'WindowButtonDownFcn', {@after_click, gobj, iniprop});
    
%     if first
%         set_prop(gobj, iniprop)
%     end
    set(gcf, 'WindowButtonUpFcn', @dispname);
end

% function first = isfirst
% ao = findobj(gcf, '-property', 'displayname');
% adisp = {ao(:).DisplayName};
% stit = get(get(gca,'title'),'string');
% if ~any(strcmp(adisp, stit))
%     first = true;
% else
%     first = false;
% end

function after_click(~, ~, gobj, prop)
set_prop(gobj, prop)

        
function set_prop(gobj, prop)
Np = length(prop(:,1));
for i = 1 : Np
    set(gobj, prop{i, 1}, prop{i, 2});
end
function plot_vol(atlas, idx)
% Just to check volume
Nlab = length(atlas.tissuelabel);
if nargin > 1 && isempty(idx)
    idx = 1 : Nlab;
else
    Ni = length(idx);
    if Ni==Nlab % Logical indexation
        Nlab = sum(idx);
        idx = find(idx==1);
    else
        Nlab = length(idx);
    end
end
%---- Define atlas coordinates as a N*3 matrix (N = prod(dim) = 902629) 
dim = atlas.dim;
[X, Y, Z]  = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
Apos   = [X(:) Y(:) Z(:)];

%---- Apply transform to put atlas position in head coordinates
% Left and Right become like fieldtrip and BS convention
Mtr = atlas.transform; 
pos = Apos;
pos(:,4) = 1;
pos = pos * Mtr';
pos = pos(: , 1:3);

%---- Associate tissue identification number
id_tis = atlas.tissue(:); %  902629x1 

%---- Figure
Alab = atlas.tissuelabel;

Ntis = length(Alab);
colc = color_group(Ntis);

figure, 
set(gcf,'units','centimeters', 'position', [10 7 28 20],'color',[0 0 0])
set(gca,'color',[0 0 0], 'position',[0.005 0.00 .99 .92])

hold on

for j = 1 : Nlab 
    ig = find(id_tis == idx(j));
    slab = Alab{j};

    plot3(pos(ig,1), pos(ig,2), pos(ig,3),'o','markersize',2,...
        'markerfacecolor',colc(j,:),'markeredgecolor','none',...
        'displayname', slab)
end
view(90, 90)
lightangle(140, 50)
axis tight equal off;
set(gcf, 'WindowButtonDownFcn', @dispname);
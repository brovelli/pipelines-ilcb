function plist = make_pathlist(Carch)
% MAKE_PATHLIST Genere un ensemble de chemins d'acces aux dossiers selon
% l'architecture definie dans la cellule archicell
% Ces chemins sont stockes dans la cellule pathlist.
%
% Ex. d'architecture :
%
%     D1              D2        D3        D4       
%
%                                   |--- Run_1
%              	|---- S01 |---- MEG |--- Run_2
%               |                   |--- Run_3
%               |
%               |---- S02 |---- MEG |--- Run_2    
%               |                   |--- Run_3
% Project_1 ----|                  
%               |                                 
%              	|---- P01 |---- MEG |--- Run_3
%               |
%               |---- P02 |---- MEG |--- Run_1    
%                                   |--- Run_2
%
% Construction de l'argument archicell :
% Une cellule comportant autant de ligne que de niveau n de dossier (Dn)
% et 2 colonnes : col. 1 : le nom complet ou partiel du dossier (le debut
%                 du nom a recherche)
%                 col. 2 : l'indication sur le caractere complet ou partiel
%                 du nom du dossier cherche (1 : nom partiel, 0 : complet)
% Ex. avec l'architecture ci-dessus, archicell = p avec :
% p{1,1}= {'Project_1'};  	p{1,2}= 0;
% p{2,1}= {'S','P'};       	p{2,2}= 1;
% p{3,1}= {'MEG'};          p{3,2}= 0;
% p{4,1}= {'Run'};          p{4,2}= 1;
% 
% p{2,2}= 1 => l'ensemble des dossiers dont le nom commence par 'S' ou 'P' 
% sera recherche a l'interieur du dossier D2 pour definir les chemins 
% d'acces. Idem avec p{4,1}, tous les dossiers dont le nom commence par
% 'Run' trouves dans le dossier D4 seront integres a la liste de chemins.
% 
% 
% Si le chemin du dossier de travail courant de Matlab n'est pas celui
% ou sont situes les dossiers en D1, il est necessaire d'indiquer le chemin
% d'acces complet au dossier place en position D1 dans p{1,1}.
% Ex. : p{1,1}= {'C:\path1\Project_1'};
%
% >> pathlist = make_pathlist(p)
% >> char(pathlist)
% >> ans =
%   F:\path1\Project_1\S01\MEG\Run_1 
%   F:\path1\Project_1\S01\MEG\Run_2
%   F:\path1\Project_1\S01\MEG\Run_3
%   F:\path1\Project_1\S02\MEG\Run_2
%   F:\path1\Project_1\S02\MEG\Run_3
%   F:\path1\Project_1\P01\MEG\Run_3
%   F:\path1\Project_1\P02\MEG\Run_1
%   F:\path1\Project_1\P02\MEG\Run_2
      
plist = [];

Carch = check_arch(Carch);
if isempty(Carch)
    return;
end

% Initialisation des chemins initiaux
inipath = cell(1000, 1);
k = 1;

Ni = length(Carch{1, 1});
for i = 1 : Ni
    [pini, dnam] = fileparts(Carch{1,1}{i});
    dlist = dirlist(pini, {dnam, Carch{1,2}});
    if ~isempty(dlist)
        Nd = length(dlist);
        for j = 1 : Nd
            inipath{k} = dlist{j};
            k = k + 1;
        end
    end
end
inipath = inipath(~cellfun(@isempty, inipath));
Na = length(Carch(:, 1));
if Na==1 && ~isempty(inipath)
    plist = inipath;
    return
end
Ns = length(inipath);
% Find folder matching with Carch(2, 1) name list
for j = 1 : Ns
    pini = inipath{j};
    dlist = dirlist(pini, Carch(2, :));
    if ~isempty(dlist)
        if j==1
            parta = dlist;
        else
            parta = [parta ; dlist]; %#ok
        end
    end
end
if ~exist('parta', 'var')
    warning('--- Directories not found with input folder structure');
    return;
end

plist = parta;

if Na > 2
    Carch = [{ parta, 0};  Carch(3:end, :)];
    plist = make_pathlist(Carch);
end

% Make list of existing directories matching with Cfold(1) cell of
% strings
function dlist = dirlist(dpath, Cfold)
opt = Cfold{2};
cdir = Cfold{1};
if ischar(cdir)
    cdir = {cdir};
end
Nd = length(cdir);

% Initialize path list
dlist = cell(1000,1);

k = 1;
for i = 1 : Nd
    dnam = cdir{i};
    if opt
        dnam = strrep([dnam, '*'], '**', '*');
        pdir = [dpath, fsep, dnam];
        ddir = dir(pdir);
        if ~isempty(ddir)
            Ns = length(ddir);
            for j = 1 : Ns
                % Be sure it's a folder
                pdirf = [dpath, fsep, ddir(k).name];
                if exist(pdirf, 'dir')
                    dlist{k} = pdirf;
                    k = k + 1;
                end
            end
        end
    else
        pdir = [dpath, fsep, dnam];
        if exist(pdir, 'dir')
            dlist{k} = pdir;
            k = k + 1;
        end
    end
end
dlist = dlist(~cellfun(@isempty, dlist));

function Carch = check_arch(Carch) 
if isempty(Carch) || (iscell(Carch) && isempty(Carch{1}))
    Carch = [];
    return;
end

Na = length(Carch(:, 1));
% Be sure searching strings are in a cell
for i = 1 : Na
    if ischar(Carch{i, 1})
        Carch{i, 1} = Carch(i, 1);
    end
end
% Remove the blanks at the beginning or end of the search strings
Carch(:, 1) = cellfun(@(x) regexprep(x, '^\s* | \s*$', ''), Carch(:, 1), 'UniformOutput', 0);
% Remove the final filesep if any
Carch(:, 1) = cellfun(@(x) regexprep(x, ['\', fsep, '$'], ''), Carch(:, 1), 'UniformOutput', 0);

% Initialisation paths
pth_ini = Carch{1,1};
Ni = length(pth_ini);
for i = 1 : Ni
    pini = pth_ini{i};
    if ~containss(pini, fsep) 
        Carch{1,1}{i}=[pwd, fsep, pini];
    end
end

      
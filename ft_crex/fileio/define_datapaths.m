function [allpaths, subj, grp, proj] = define_datapaths(pcell, isubj, igrp, iproj)
% In addition to return all the data paths found using pcell directories
% architecture (see make_pathlist), the subject and group names are
% associated to each data path 
% Works only if the subject and group names appear as directory name in
% directories architecture.
% isubj : position in pcell of the directory corresponding with subject name 
% igrp : position of the directories holding group name according to pcell 
% iproj : position of the directory layer in pcell representing project's
% name
% 
% Example: data are stored in BaPa project directory
% It is split in 2 groupes of participants CAC and DYS
% Inside each group directory, the subject directories are nammed after the
% subject name (S01, S02, S03...). In this configuration:
%
% pcell = {{'F:\bapa'} , 0
%           {'CAC', 'DYS'}, 0 
%           {'S'}, 1
%          }; 
% with isubj = 3, igrp = 2 and iproj = 1
%

% Define all path list
allpaths = make_pathlist(pcell);

if isempty(allpaths)
    [subj, grp, proj] = deal([]);
    return;
end

Np = length(allpaths);

% Use partpath function to find the corresponding subject, group and
% project names
if nargin>=2 && ~isempty(isubj) 
    subj = partpath(pcell, allpaths, isubj);
else
    subj = cell(Np, 1);
end

if nargin>=3 && ~isempty(igrp)
    grp = partpath(pcell, allpaths, igrp);
else
    grp = cell(Np, 1);
end

if nargin==4 && ~isempty(iproj)
    proj = partpath(pcell, allpaths, iproj);
else
    proj = cell(Np, 1);
end
                
function partlist = partpath(pcell, allpaths, ipart)    

Nall = length(allpaths);    
partlist = cell(Nall, 1);

if ~isempty(ipart)
    
    charall = char(allpaths);

    partpath = make_pathlist(pcell(1:ipart,:));

    Npart = length(partpath);

    for i = 1 : Npart
        pattern = partpath{i};
        [~, name] = fileparts(pattern);
        ind = strcmp(pattern, cellstr(charall(:, 1:length(pattern))));
        
        partlist(ind) = {name};
    end
end
        
    


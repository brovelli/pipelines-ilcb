function Meach = read_trans_each(pref, subjname, pdb_bv, pdb_fs)
% Read the successive transformation matrix files from BrainVisa processing 
% that are indicated in the referential.txt file
% Input :
% - pref : referential TXT file path [REQUIRED]
% - subjname : subject's name [REQUIRED]
%
% If the database is no longer at the same location as when processed by Brainvisa 
% and Freesurfer, it is possible to indicate the new path to the databases:
% - pdb_bv : path to the brainvisa database [default : [] (use the path in the
% referential txt file]
% - pdb_fs : path to the freesurfer database [default: []]
%
% Output : the total affine transform matrix
%
% TO DO: allow more flexibility in database name
% => chercher le dernier dossier commun entre le nom du chemin en input et le
% nom dans le path dans referential.txt

db_bv = 'db_brainvisa';
db_fs = 'db_freesurfer';

fid = fopen(pref, 'r');
if fid < 1
    fprintf('\n!!!\n Expected referential.txt path: %s\n', pref);
    warning('referential.txt not found');
    Meach = [];
    return;
end
ctrans = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

ctrans = ctrans{1};
Nt = length(ctrans);

% 
% Replace {0} by the subject's name
ctrans = strrep(ctrans, '{0}', subjname);

Mtot = eye(4, 4);
vs = [0 0 0 1];
Meach = cell(Nt, 1);

% Read each transform matrix
for i = 1 : Nt
    strans = ctrans{i};
    % Use the path to subject's brainvisa or freesurfer database provided in
    % input whereas the path in referential.txt file
    isinv = strcmp(strans(1:3), 'inv');
    
    ibv = strfind(strans, db_bv);
    ifs = strfind(strans, db_fs);
    if ~isempty(ibv)
        % fullfile will change the filesep direction if the OS for current data
        % processing is different from the OS used for atlas preparation
        spath = fullfile(pdb_bv, strans(ibv(end) + length(db_bv) : end));
    elseif ~isempty(ifs)
        spath = fullfile(pdb_fs, strans(ifs(end) + length(db_fs) : end));       
    end
    if ~exist(spath, 'file')
        warning('Unable to access to the transform matrix file\n%s\n', spath);
        Meach = [];
        return;
    else
        fid = fopen(spath, 'r');
        Cmat = textscan(fid, '%f');
        fclose(fid);
        Mf = reshape(Cmat{1}, 3, 4)';
        Mt = [[Mf(2:end, :) Mf(1, :)'] ; vs];
        if isinv
            Mtot = Mtot / Mt;
            Meach{i} = inv(Mt);
        else
            Mtot = Mtot * Mt;
            Meach{i} = Mt;
        end
        
    end
    
end




        



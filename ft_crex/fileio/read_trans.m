function [Meach, Mtot] = read_trans(opt)
% Read the successive transformation matrix files from BrainVisa/Freesurfer processing 
% as indicated in the referential.txt file
% 
% --- Inuput: opt structure with fields:
%
%   - ref: path to the referential txt file [default: the one inside ft_crex
%       toolbox: ft_crex/atlas/referential.txt]
%
%   - trans_dir: path of a directory containing all the required trm files
%       [default: [] - trm are read directly from the brainvisa and freesurfer
%       databases according to opt.db_dir paths
%
%   - db_dir: path to the subject's brainvisa and freesurfer databases (db_dir.bv and db_dir.fs)
%       to find the trm files as indicated in opt.ref file
%       If opt.trans_dir is set, theses database paths are not required as trm
%       files have been previously imported (cf. by cp_db_import)
%
%
% --- Output: 
%
%   - Meach: cell holding each affine transform matrix extracted from trm files according to opt.ref file
%   
%   - Mtot: the total transformation matrix from the successive trm
%       transformation according to opt.ref file
%
%-CREx-180702

dopt = struct('ref', fullfile(ptool, 'atlas', 'referential.txt'), ...
            'trans_dir', '',...
            'db_dir', []);
        
opt = check_opt(opt, dopt);


fid = fopen(opt.ref, 'r');
if fid < 1
    fprintf('\n!!!\n Expected referential.txt path: %s\n', opt.ref);
    warning('referential.txt not found');
    Meach = [];
    return;
end

ctrans = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

% Determine if TRM files have been already imported inside a opt.trans_dir
% folder
if ~isempty(opt.trans_dir) && ~isempty(dir([opt.trans_dir, filesep, '*.trm']))
    isimp = 1;
    ptrans = opt.trans_dir;
else
    % If not, try to find the TRM files from the BRAINVISA and FREESURFER
    % databases according to referential.txt file
    isimp = 0;
end

ctrans = ctrans{1};
Nt = length(ctrans);

Mtot = eye(4, 4);
vs = [0 0 0 1];
Meach = cell(Nt, 1);

% Read each transform matrix
for i = 1 : Nt
    strans = ctrans{i};
    % Use the path to subject's brainvisa or freesurfer database provided in
    % input whereas the path in referential.txt file
    isinv = strcmp(strans(1:3), 'inv');
    
    if ~isimp
        isi = strfind(strans, '/{0}/');
        if ~isempty(strfind(strans, 'brainvisa'))  %#ok
            % fullfile will change the filesep direction if the OS for current data
            % processing is different from the OS used for atlas preparation
            ptr = fullfile(opt.dir.bv, strans(isi+5 : end));
        else
            ptr = fullfile(opt.dir.fs, strans(isi+5 : end));    
        end
        
    else
        % trm file should be find inside opt.trans_dir folder
        [~, ftrm, ext] = fileparts(strans);
        ptr = fullfile(ptrans, [ftrm, ext]);
    end
    ptr = strrep(ptr, '{0}', '*');
    dtr = dir(ptr);
            
    if isempty(dtr)
        warning('Unable to access the transform matrix file\n%s\n', ptr);
        Meach = [];
        return;
    end
    ptr = [dtr.folder, filesep, dtr(1).name];     
        
    fid = fopen(ptr, 'r');
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
 
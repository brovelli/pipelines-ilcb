function filtData = cmeg_preproc_filt(dpath, opt)
% Apply filter on continuous dataset
% 
% dpath: path of directory where dataset is stored (mat-file from previous processing)
%
% Parameters of the filtering option in opt structure:
%
% opt.type : kind of filter identified by a string
%       - 'hp' : high-pass filter
%       - 'lp' : low-pass filter
%       - 'bp' : band-pass filter
%       - 'ask' : ask for it, for each dataset [default]
%
% opt.fc : cut-off frequency vector in Hz
%       - size 1x2 for band-pass [f_low f_high] (ex. : [0.5 300])
%       - size 1x1 for low or high-pass filter 
%           (ex. filtopt.type = 'lp'; filtopt.fc = 300)
%       [default : [] - no filtering]
%
% opt.datatyp : prefix of mat-dataset to load and filter in data directory
%       - 'raw' : search for "rawData*.mat" data files in path directory
%       - 'filt' : search for data previously filtered "filtData*.mat, and
%                   re-apply filter on it
%       - 'clean' : dataset of data already cleaned by ICA components
%       rejection ("cleanData*.mat")
%       - any custom prefix string (search for "[custom]Data*.mat" on path
%       directory)
%       [default : 'raw']
%
% figopt.figflag : draw the first 3 channels continuous data in subplot
%           of before and after filtering to see the global effect
%           1 : draw and save [default] ; 0 : no figure
%
% Save the new filtering data in dpath directory
%
%____
%-CREx 20140520 

fprintf('\n\t\t-------\nApply filter on dataset\n\t\t-------\n')
fprintf('\nProcessing of data in :\n%s\n\n', dpath);

dopt = struct('type', 'ask',...
            'fc', [],...
            'figflag', 1,...
            'savepath', dpath,...
            'figpath', [],...
            'channel', '*');
        
dopt.datatyp = {'raw', 'clean', 'filt'};

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% Check for data
[pmat, nmat] = find_mat(dpath, opt);

% Check for filt options
opt = cmeg_filt_opt(opt);
if (isempty(opt.type) || strcmp(opt.type, 'none')) 
    if ischar(opt.channel)
        filtData = [];
        return;
    end
end

% Try to extract data from raw 4D or EEG file
if isempty(pmat)
    rawData = cmeg_extract_raw(dpath);
    pmat = 'RAW-file';
    nmat = 'rawData';
    vnam = nmat;
    if isempty(rawData)
        warning('No data found')
        filtData = [];
        return;
    end
else    
    [rawData, vnam] = loadvar(pmat);
end


fprintf('\nInput data :\n%s\n', pmat);
    
% Apply filter
filtData = cmeg_filt(rawData, opt);

if ~nargout
    % Get the preproc suffix indication
    strp = preproc_str(opt);

    if ~isfield(filtData.hdr, 'info')
        filtData.hdr.info.type = 'continuous_filt';
    else
        styp = filtData.hdr.info.type;
        filtData.hdr.info.type = strrep(styp, 'raw', 'filt');
    end

    % Keep preprocessing info
    filtData = add_info(filtData, 'sfilt', strp);

    %- Save data in dpath directory if nargout==0
    nmat = strrep(nmat, 'raw', 'filt');
    savnam = cmeg_name(nmat, strp);

    % Save path
    spath = opt.savepath;
    pfilt = [spath, filesep, savnam];
    
    % Keep initial data structure name but change 'raw' to 'filt'
    snam = strrep(vnam, 'raw', 'filt');
    
    S = [];
    S.(snam) = filtData;    
    
    save(pfilt, '-v7.3', '-struct', 'S');
    fprintf('\nFiltered data set saved as: %s\n', savnam);
end

if opt.figflag
    if isempty(opt.figpath)
        pfig = make_dir([spath, filesep, 'filt', strp], 1);
    end
    
    % Make some figures of the filtered vs raw data
    opt.savepath = pfig;
    opt.chan = filtData.label(1:3);
    opt.info = pfilt;
    cmeg_filt_fig(filtData, rawData, opt) 
end


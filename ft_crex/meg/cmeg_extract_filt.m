function filtData = cmeg_extract_filt(dpath, opt)
% Extract and apply filter to continuous data set 
% 
% dpath: path of directory of the raw data set
%
% Parameters of the filtering option in opt structure:
%
% opt.type : kind of filter identified by a string
%       - 'hp' : high-pass filter
%       - 'lp' : low-pass filter
%       - 'bp' : band-pass filter
%       - 'ask' : ask for it, for each dataset
%       - 'none' : no filter [default]
%
% opt.fc : cut-off frequency vector in Hz
%       - size 1x2 for band-pass [f_low f_high] (ex. : [0.5 300])
%       - size 1x1 for low or high-pass filter 
%           (ex. filtopt.type = 'lp'; filtopt.fc = 300)
%       [default : [] - no filtering]
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

dopt = struct('type', 'none',...
            'fc', [],...
            'figflag', 1,...
            'savepath', dpath,...
            'figpath', [],...
            'channel', '*');

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% Check for filt options
opt = cmeg_filt_opt(opt);
if isempty(opt.type) || strcmp(opt.type, 'none')
    israw = 1;
else
    israw = 0;
end

% Extract data from raw 4D or EEG file
rawData = cmeg_extract_raw(dpath);

if isempty(rawData)
    warning('No data found')
    filtData = [];
    return;
end

% Apply filter / remove bad channels even if filter is not set
filtData = cmeg_filt(rawData, opt);

% Automatically correct the 1s at the edge of the continuous data to remove
% filter artefact (that appears even with a 10s border padding method)
% No stimulation is expected so close to the border
tend = filtData.time{1}(end);
filtData = cmeg_artefact_rm(filtData, [0 2; tend-2 tend]);

if ~nargout
    if israw
        snam = 'rawData';
    else
        snam = 'filtData';
    end
    % Get the preproc suffix indication
    strp = preproc_str(opt);

    filtData.hdr.info.type = 'continuous_filt';

    if ~israw
        % Keep preprocessing info
        filtData = add_info(filtData, 'sfilt', strp);

        %- Save data in dpath directory if nargout==0
        savnam = cmeg_name(snam, strp);
    else
        savnam = snam;
    end

    % Save path
    spath = opt.savepath;
    pmat = [spath, filesep, savnam];
    
    S = [];
    S.(snam) = filtData;    
    
    save(pmat, '-v7.3', '-struct', 'S');
    fprintf('\nFiltered data set saved: %s\n', pmat);
end

if opt.figflag && ~israw
    if isempty(opt.figpath)
        pfig = make_dir([opt.savepath, filesep, 'filt', preproc_str(opt)], 1);
    end
    
    % Make some figures of the filtered vs raw data
    opt.savepath = pfig;
    opt.chan = filtData.label(1:3);
    opt.info = dpath;
    cmeg_filt_fig(filtData, rawData, opt) 
end


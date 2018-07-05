function hs = read_headshape(pmeg)
% Try to find heashape file in pmeg directory 
% for 4D and Neuromag system



phs = fullfile(pmeg, 'hs_file');

if exist(phs, 'file')
    hs = ft_read_headshape(phs, 'unit', 'mm');
    return;
end

% Try to find raw MEG data (cf. for fif files that hold headshape)
hs = [];
praw = filepath_raw(pmeg);
if ~isempty(praw)
    try
        hs = ft_read_headshape(praw, 'unit', 'mm');
    catch
        fprintf('No headshape file found in %s directory', pmeg)
    end
end
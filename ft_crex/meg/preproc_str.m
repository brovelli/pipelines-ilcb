function strproc = preproc_str(popt)
% Expecting preprocessing suffix string depending on preprocessing option
% done on epoched data
%
%

strproc = '';
if isempty(popt)
    return;
end


%---- FILTER 
if isfield(popt, 'type') && isfield(popt, 'fc') && ~isempty(popt.fc)
    fc = popt.fc;    
    % HIGH-PASS OR LOW-PASS
    if strcmp(popt.type, 'hp') || strcmp(popt.type, 'lp')
        fcs = num2str(fc);
        fcs(fcs=='.') = '';
        strproc = [strproc,'_',popt.type, fcs];
    % BAND-PASS
    elseif strcmp(popt.type, 'bp')
        fcs = ['hp', num2str(fc(1)),'_lp', num2str(fc(2))];
        fcs(fcs=='.') = '';
        strproc = [strproc,'_', fcs];
    end
end

%---- RESAMPLE
if isfield(popt, 'res_fs') && ~isempty(popt.res_fs) && popt.res_fs > 0   
    fss = num2str(popt.res_fs);
    fss(fss=='.') = '';
    strproc = [strproc,'_res',fss];
end
%---- CROP WINDOW
if isfield(popt, 'crop_win') && ~isempty(popt.crop_win) && sum(abs(popt.crop_win)) > 0   
    dur = num2str(sum(abs(popt.crop_win)));
    if dur(1)=='0'
        dur(dur=='.') = '';
    else
        dur(dur=='.') = 'p';
    end       
    strproc = [strproc,'_crop',dur,'s'];
end  

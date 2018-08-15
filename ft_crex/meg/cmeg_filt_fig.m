function cmeg_filt_fig(filtData, rawData, opt)
% Figures of subplots showing filtered data vs raw data (one figure / channel)
%
% -- Inputs:
% * filtData : fiedtrip data structure of filtered data
% * rawData : fieldtrip data structure of raw data
%
% * opt : initial parameters structure for filtering
%   opt.type: filter type ('lp', 'hp', 'bp')
%   opt.fc : associated cut-off frequency
%   
%   opt.chan: list of the channels of the data to plot [default: all in filtData.label]
%   opt.savepath: directory where figures will be saved [default: pwd, filesep, 'filt_fig']
%   opt.info: additionnal information to put on the figure title (ex. data path)
%
%____
%-CREx-180427

% Check for input option
dopt = struct('savepath', [], 'type', '', 'fc', [], 'chan', [], 'info', '');

if nargin < 3 
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

pfig = opt.savepath;
if isempty(pfig) || ~exist(pfig, 'dir')
    pfig = make_dir([pwd, filesep, 'filt_fig'], 1);
end

%-- Define labels specifying kind of filter
[finfo, ssav] = filt_str(opt);

chan = opt.chan;
if isempty(chan)
    chan = filtData.label;
end

% Number of data to plot
Nc = length(chan);

% Try to define data units
sun = det_units(chan{1});

%-- Make figure for each channel
for i = 1 : Nc
    slab = chan{i};
    
    % Get time and amplitude vectors 
    [tfilt, xfilt] = extract_data(filtData, chan{i});    
    [traw, xraw] = extract_data(rawData, chan{i});
    
    % Title string
    stit = {['Filtering result: ', chan{i}]; opt.info};
    
    figure
    set(gcf, 'visible', 'off', 'units','centimeter','position',[1 1 27 17])
    
    %-- Filtered data subplot
    subplot(211)
    plot_dat(tfilt, xfilt, sun);
    
    % Add title
    title(stit, 'fontsize',14, 'interpreter', 'none')
    
    %-- Raw data subplot
    subplot(212)
    plot_dat(traw, xraw, sun); 
    
    % Add filtering indication + move subplot
    format_plot(finfo)
       
    % Save the figure
    % Name is based on the frequency cut-off parameters and the channel 
    export_fig([pfig, filesep, 'filt', ssav,'_', slab, '.png'], '-m1.5')
    
    close
end

% Get filtering indications as string to add to the figure (figsav) as well
% as suffix to add to the jpeg file that is saved (savstr)
function [finfo, ssav] = filt_str(opt)
% Filter characteristics according to options structure
% filter type ('bp','lp' or 'hp')
ftyp = opt.type;
% Frequency cut-off in Hz
fcut = opt.fc;
if ~isempty(ftyp)
    if numel(fcut) == 1
        sfc = num2str(fcut);
        finfo = ['High-pass: f_{c} = ', sfc,' Hz'];
        if strcmpi(ftyp, 'lp')
            finfo = strrep(finfo, 'High', 'Low');
        end
    else
        sf1 = num2str(fcut(1));
        sf2 =  num2str(fcut(2));
        finfo = ['Band-Pass: [', sf1,' - ', sf2,' ] Hz'];
    end
end
ssav = preproc_str(opt);
if isempty(ssav)
    ssav = 'filt';
end

function sun = det_units(lab)
if strcmp(lab(1), 'A') || strcmpi(lab(1), 'M')
    sun = 'T';
else
    sun = 'V';
end

function format_plot(finfo)
% Format labels if needed (add x10^n indication in label string)
format_label

% Add indication on the figure
put_figtext('B) Original data', 'nw', 12, 'w', 'o');
% Move the subplot down
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2)+0.05 pos(3:4)])

subplot(211)
put_figtext(['A) ', finfo], 'nw',12, 'w', 'o');
format_label 

function plot_dat(tdat, xdat, sun)
plot(tdat, xdat);
xlim(tdat([1 end]))

% Add axes label and change fontsize
xlabel('Time (s)','fontsize',14)
ylabel(['Amplitude (', sun, ')'],'fontsize',14)

set(gca,'fontsize',14)

% Extract time and amplitude vector of the request channel (chanstr)
function [tdat, xdat] = extract_data(ftdat, chanstr)
tdat = ftdat.time{1};
xdat = ftdat.trial{1}(strcmp(ftdat.label, chanstr)==1, :);


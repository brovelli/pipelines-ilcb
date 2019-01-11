function cp_inv_dynmesh_fig(Spow, opt)

% Default option
dopt = struct('decim', 1,...
            'xlim', Spow.time([1 end]),...
            'savepath', pwd,...
            'clim', [],...
            'info', '',...
            'title', '',...
            'output', 'both',...
            'condition', '');
% output: 'avi', 'gif' or 'both'
if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end
 
if ~isempty(opt.info)
    opt.info = ['_', opt.info];
end

dc = opt.decim;
xl = opt.xlim;
time = Spow.time;
if ~isempty(xl)
    ind = find(time >= xl(1) & time <= xl(2));
else
    ind = 1 : length(time);
end
ind = ind(1) : dc : ind(end);

Nt = length(ind);

pfig = opt.savepath;
% Keep png frame
pfram = make_dir([pfig, fsep, 'frames']);

% Output video or GIF  
otyp = opt.output;
if strcmp(otyp, 'both') || strcmp(otyp, 'gif')    
    gfig = [pfig, fsep, 'powcomp', opt.info, '.gif'];
    isgif = 1;
else
    isgif = 0;
end
dt_frame = diff(time(ind(1:2)))*16;

if strcmp(otyp, 'both') || strcmp(otyp, 'avi')
    vid = VideoWriter([pfig, fsep, 'powcomp', opt.info, '.avi']);
    vid.FrameRate = round(1/dt_frame);
    vid.Quality = 100;
    open(vid);
    isvid = 1;
else
    isvid = 0;
end

pow = Spow.pow;
Smesh = Spow.mesh;

if isempty(opt.clim)
    mp = [min(pow(:)) max(pow(:))];
    dc = diff(mp);
    opt.clim = [mp(1)+0.2*dc mp(2)-0.3*dc];
end
wb = waitbar(0, 'Video of cortical activity...', 'name', 'Dynmesh video');
wb_custcol(wb, [0 0.05 0.93]);
for j = 1 : Nt
    
    waitbar((j-1)/Nt, wb, ['Video:', opt.info, ' - frame: ', num2str(j)]);
    
    ii = ind(j);
    stim = num2str(time(ii)*1e3, '%3.0f');
    st = ['t = ', stim, ' ms'];
    opt.timer = [st blanks(length('t = -1000.0 ms') - length(st))]; 

    powt = pow(:, ii);

    hfig = dynmesh_plot(Smesh, powt, opt);

    drawnow
    % Capture the plot as an image 
    frame = getframe(hfig);  
    [imi, cm] = rgb2ind(frame2im(frame), 256); 

    % Write to the GIF File 
    if isgif
        if j == 1
            imwrite(imi, cm, gfig, 'gif', 'Loopcount', Inf, 'Delaytime', dt_frame);
        else
            imwrite(imi, cm, gfig, 'gif','WriteMode', 'append', 'Delaytime', dt_frame);
        end
    end
    
    if isvid
        writeVideo(vid, frame)
    end
    
    % pause(0.001)
    snum = num2str(j);
    sfram = [repmat('0', 1, 3 - length(snum)) snum];
    fnam = ['powcomp', opt.info, '_', sfram, '_', name_save(stim),'ms'];
    export_fig([pfram, filesep, fnam, '.png'], '-m1.2')
    close(hfig)
end
close(wb);
if isvid
    close(vid);
end

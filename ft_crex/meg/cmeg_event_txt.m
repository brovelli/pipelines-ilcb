function cmeg_event_txt(hdr_event, ptxt, praw)
% Create a txt file with data information and list of events (markers)
%-CREx181120 

%________
% Check input
if nargin < 3 
    praw = '';
end

% ________
% Prepare data to display

Sev = hdr_event.event;
fevent = fieldnames(Sev);

% All events: cell of char for print
S = struct2cell(Sev)';
Ss = cellfun(@num2str, S, 'UniformOutput', 0);
Nr = length(Ss(:, 1));
Sc = char(Ss);
nch = length(Sc(1, :));
Ss = cellfun(@(x) [x blanks(nch-length(x))], Ss, 'UniformOutput', 0);
fev = cellfun(@(x) [x blanks(nch-length(x))], fevent, 'UniformOutput', 0);

% Additionnal info in header of ptxt file
hdr = hdr_event;
hdrd = struct('nChans', NaN, 'nSamples', NaN, 'Fs', NaN);
hdr = check_opt(hdr, hdrd);

% Count the diffenrent type of events
Styp = [];
if any(strcmp(fevent, 'type')) && any(strcmp(fevent, 'value'))
    [utyp, ~, iu] = unique({Sev(:).type}');
    Nu = length(utyp);
    for i = 1 : Nu
        ftyp = utyp{i};
        vals = {Sev(iu==i).value}';
        vals = cellfun(@num2str, vals, 'UniformOutput', 0);
        
        [uvals, ~, uv] = unique(vals);
        Nv = length(uvals);
        Styp.(ftyp).strval = uvals;
        nval = zeros(Nv, 1);
        for j = 1 : Nv
            nval(j) = sum(uv==j);
        end
        Styp.(ftyp).numval = nval;
    end
end

% ________
% Display data informations & events 
fid = fopen(ptxt, 'w');

if ~fid
    warning('Unable to create event-info txt file as %s', ptxt)
    return;
end

fprintf(fid, '\n\t\t--------\n\t\tList of events\n\t\t--------\n\n');
fprintf(fid, '\nData path: %s\n\n', praw);   
fprintf(fid, 'Number of channels: %d\n', hdr.nChans);
fprintf(fid, 'Recording duration (s): %4.1f\n', (hdr.nSamples - 1)./hdr.Fs);
fprintf(fid, 'Sample frequency (Hz): %4.1f\n', hdr.Fs);
fprintf(fid, '\n       --------------------------\n');
fprintf(fid, ['\n\t\t', strjoint(fev, '\t\t'), '\n\n']);
for i = 1 : Nr
    fprintf(fid, ['\t\t', strjoint(Ss(i, :), '\t\t'), '\n']);
end
fprintf(fid, '\n\n      --------------------------\n\n');

if ~isempty(Styp)
    fprintf(fid, '\n %d type(s) of event found:\n', Nu);
    for i = 1 : Nu
        ftyp = utyp{i};
        fprintf(fid, '\n-- Type n°%d: %s with values:\n', i, ftyp);
        fprintf(fid, '\n\t\tValues\t(Number)\n\n');
        vals = Styp.(ftyp).strval;
        Nv = length(vals);
        for j = 1 : Nv
            fprintf(fid, '\t\t%s\t(%d)\n', vals{j}, Styp.(ftyp).numval(j));
        end
    end
end
fprintf(fid, '\n      --------------------------\n');
fclose(fid);

        
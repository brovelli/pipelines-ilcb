function opt = cmeg_filt_opt(opt)
% Check for the input filt option opt structure
% Ask for filter type and frequency cut-off if opt.type=='ask' or if frequency
% cut-off vector (opt.fc) is empty or incorrectly set
%

dopt = struct('type', 'ask', 'fc', [], 'fig', 1);
if nargin < 1 || isempty(opt)
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

fcin = opt.fc;
if ~isempty(fcin)
    fcin = fcin(fcin > 0);
end

if (~strcmp(opt.type, 'ask') && isempty(fcin)) ||...
        ((strcmp(opt.type, 'lp') || strcmp(opt.type, 'hp')) && length(fcin)~=1) ||...
        (strcmp(opt.type, 'bp') && length(fcin) ~= 2)
    opt.type = 'ask';
end
    
% Ask mode: filter type and cut-off frequencies are asked at the command prompt
if strcmp(opt.type,'ask')
    fprintf('\n\t\t--------\n\tFiltering options\n\t\t--------\n\n');
    disp('Choose the filter to apply:')
    disp('   None       -> 0')
    disp('   High-pass  -> 1')
    disp('   Low-pass   -> 2')
    disp('   Band-pass  -> 3')
    onum = input('              -> ');
    if ~onum
        opt.type = 'none';
        return
    end
    fcz = zeros(1,2);
    
    disp(' ')
    if onum==1 || onum==3
        fcz(1) = input('High-pass cut-off frequency (Hz): ');
        if onum==1
            opt.type = 'hp';
        else
            opt.type = 'bp';
        end
    end
    
    if onum==2 || onum==3
        fcz(2) = input('Low-pass  cut-off frequency (Hz): ');
        if onum==2
            opt.type = 'lp';
        end
    end
    disp(' ')
    opt.fc = fcz(fcz > 0);
end




            

            
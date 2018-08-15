function fopt = check_ica_filt_opt(fopt, fc_hp)
%% TO DO: see if it is really required to have the same HP as those used for ICA on 
% the resampled data version for ICA cleaning
if isempty(fopt.type) || strcmp(fopt.type, 'lp')
    fprintf('\n-----\n');
    warning('for ICA rejection a High-Pass filter is required to remove the LF components')      

    if isempty(fopt.type) 
        fopt.type = 'hp';
        fopt.fc = fc_hp;
        fprintf('\n --> HP filter with fc = %1.1f Hz will be applied\n\n', fc_hp);
    elseif strcmp(fopt.type, 'lp')
        fopt.type = 'bp';
        fopt.fc = [fc_hp fopt.fc];   
        fprintf('\n --> BP filter with fc = [%1.1f %3.1f] Hz will be applied\n\n', fc_hp, fopt.fc);
    elseif strcmp(fopt.type, 'bp')
		if fopt.fc(1) < fc_hp
			fopt.fc(1) = fc_hp;
		end
		fprintf('\n --> BP filter with fc = [%1.1f %3.1f] Hz will be applied\n\n', fopt.fc(1), fopt.fc(2));
	end
    fprintf('\n-----\n');
end
function Sdo = init_do(Sdo, ir)
% Initialize the 'do' structure to set all processing for MEG data to 0 for
% the run n° ir 
% (case when a run has been declare as bad during the MEG prep_pipeline (see
% cp_meg_prep)
%

[fdo, Nd] = get_names(Sdo);
for i = 1 : Nd
    fn = fdo{i};
    vdo = Sdo.(fn);
    vdo(ir) = 0;
    Sdo.(fn) = vdo;
end


function opt = check_opt(opt, defopt)
% Check for input options and assign default value if option missing or empty
% compare to default option defopt
%
%-- CREx2017

if isempty(opt)
    opt = defopt;
    return;
end

fopt = fieldnames(defopt);

for i = 1 : length(fopt)
    fnam = fopt{i};
    if ~isfield(opt, fnam) || isempty(opt.(fnam))
        opt.(fnam) = defopt.(fnam);
    else
        % Check for substructure parameters and assign default if isempty 
        if isstruct(defopt.(fnam))
            opt.(fnam) = check_opt(opt.(fnam), defopt.(fnam));
        end
    end
end
 
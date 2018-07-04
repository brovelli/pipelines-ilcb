function Sren = rename_struct(Sini, Cstr_ini, Cstr_end)
% Rename field names of the multi-levels structure Sini
% For any k (for k = 1 : length(Cstr_ini)), replace all occurences of 
% Cstr_ini{k} by its corresponding string, Cstr_fin{k}.
% Exemple :
% Sini is a structure with these different fields :
% Sini.CAC.Morpho.time
%                .avgROI
%                .Morpho_Seman
%         .Ortho.time
%               .avgROI
%               .Ortho_Seman
% Sini.DYS.Morpho.time
%                .avgROI
%                .Seman
%          ...
% The occurrences to replace are defined in this cellule of string :
%   Cstr_ini = {'Morpho','Ortho','Seman'};
% The corresponding substitution string  :
%   Cstr_fin = {'Morphological', 'Orthographic', 'Semantic'}; 
% Sren = rename_struct(Sini, Cstr_ini, Cstr_fin) changes the field names as :
% Sren.CAC.Morphological.time
%                       .avgROI
%                       .Morphological_Semantic
%         .Orthographic.time
%                      .avgROI
%                      .Orthographic_Semantic
% Sren.DYS.Morphological.time
%                       .avgROI
%                       .Semantic
%          ...

if ~isstruct(Sini)
    Sren = Sini;
    return;
end

if ischar(Cstr_ini)
    Cstr_ini = {Cstr_ini};
end

if ischar(Cstr_end)
    Cstr_end = {Cstr_end};
end

Nc = length(Cstr_ini);
Ns = length(Sini);

for i = 1 : Ns
    Str = Sini(i);
    
    fnames = fieldnames(Str);
    Nf = length(fnames);
    
    % Could be unique condition name or combinaison (structure of results
    % from statistical analysis for example)
    % Replace any occurence of the fieldnames at this level
    for j = 1 : Nf
        fnam = fnames{j};
        for k = 1 : Nc
            stri = Cstr_ini{k};
            strf = Cstr_end{k};
            
            iocc = find_strocc(fnam, stri, strf);    
            if isempty(iocc)
                break;
            end
            
%             if length(stri) < length(strf)
%                 newname = regexprep(fnam, ['(?!\',strf,')\',stri], strf);
%             else
%                 newname = regexprep(fnam, stri, strf);
%             end
            newname = strrep(fnam, stri, strf);
            % The name has been change
            if ~strcmp(fnam, newname)
                if sum(strcmp(newname, fnames))==0
                    Sini.(newname) = Sini.(fnam);
                    Sini = rmfield(Sini, fnam);
                    fprintf('Field name %s replace by %s\n', fnam, newname);
                    fnam = newname;
                else
                    fprintf('New field %s already exists\n', newname);
                end
            end
          
        end
        Sini.(fnam) = rename_struct(Sini.(fnam), Cstr_ini, Cstr_end);
    end
 
end
Sren = Sini;
 
function iocc = find_strocc(fname, strini, strfin)        
    iocc = strfind(fname, strini); 
    if ~isempty(iocc)
        iocc_fin = strfind(fname, strfin);
        if ~isempty(iocc_fin)
            % Shouldn't be the same occurrence as Cstr_fin{k}
            [~, ia] = setxor(iocc, iocc_fin);
            iocc = iocc(ia);
        end
    end
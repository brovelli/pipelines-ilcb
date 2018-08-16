function namsav = name_save(name)
% Format filename to be saved by deleting or replacing special characters, including points,
% blanks, dashes, or "[", "]", "/", "\"
%
% If the character "-" is found before a digit, it is assumed to correspond to
% the minus sign and will then be replaced by a "m".
% If character "." or "," is found between two digits, it is assumed to correspond to
% the point of a float number and will be then be replaced by a "p"
% Blanks are removing
% name : name to format (without the extension)

namsav = name;
namsav = regexprep(namsav, '-(?=\d)', 'm');
namsav = regexprep(namsav, '(?<=\d)\.(?=\d)', 'p');
namsav = regexprep(namsav, '(?<=\d)\,(?=\d)', 'p');
namsav = regexprep(namsav, '[/|\\\[\]\s\.,-]', '_');
while ~strcmp(namsav, strrep(namsav, '__', '_'))
    namsav = strrep(namsav, '__', '_');
end


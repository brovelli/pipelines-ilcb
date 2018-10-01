function stxt = prep_tex(stxt)
% Prepare string for latex interpreter : change '\' to '/' and add a "\" before the "_" 
% (otherwise latex makes the character that follows a subscribe)

stxt = strrep(stxt, '\', '/'); 
stxt = strrep(stxt, '_', '\_');
function save_fig(pfig)
% Save matlab figure as FIG file - be sure figure visible mode is set to 'on'
%-CREx-20180815

% Add .fig extension if not set
[~, ~, ext] = fileparts(pfig);
if isempty(ext)
	pfig = [pfig, '.fig'];
end

isv = strcmpi(get(gcf, 'visible'), 'on');
% Be sure visible mode if 'on' (to view figure when reopen it)
set(gcf, 'visible', 'on')
saveas(gcf, pfig)
if ~isv
	set(gcf, 'visible', 'off');
end
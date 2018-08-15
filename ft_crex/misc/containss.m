function isc = containss(str, pattern)
% Polyfill for Matlab < 2017
isc = ~isempty(strfind(str, pattern)); %#ok 
function Csp = strsplitt(Str, sspl)

Csp = strsplit(Str, sspl);
isc = cellfun(@(x) ~isempty(x), Csp);
Csp = Csp(isc);

    
function [fnames, Nf] = get_names(Sdat)
fnames = fieldnames(Sdat);
Nf = length(fnames);
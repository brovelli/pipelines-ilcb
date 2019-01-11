function imeg = is_meg(Sdb)
% Return an array indicating which element of Sdb has a meg field which is not
% empty == the subject with available MEG data (= rows of imeg with value == 1)
if isempty(Sdb)
    imeg = 0;
    return;
end
Ns = length(Sdb);
imeg = zeros(Ns, 1);
for i = 1 : Ns
    if ~isempty(Sdb(i).meg)
        imeg(i) = 1;
    end
end
% ! logical indexing doesn't work for structures !
imeg = find(imeg==1);
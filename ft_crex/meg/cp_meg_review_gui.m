function isval = cp_meg_review_gui(Sdb, opt)

Srev = prepare_list(Sdb);
Srev.Sdb = Sdb;
Srev.opt = opt;
ischg = preproc_review(Srev);
isval = ~ischg;
function Srev = prepare_list(Sdb)
Ns = length(Sdb);
% Id for run number 
cidr = cell(Ns, 1);
% Id for subject number
cids = cell(Ns, 1);
for i = 1 : Ns
    Nr = length(Sdb(i).meg.rundir);
    cidr{i} = (1: Nr)';
    cids{i} = repmat(i, Nr, 1);
end

isubj = cell2mat(cids);
irun = cell2mat(cidr);

Nd = length(cell2mat(cids));
dinfo = cell(Nd, 1);
for i = 1 : Nd
    Sdbs = Sdb(isubj(i));   
    dinfo{i} = [Sdbs.sinfo, '-', Sdbs.meg.rundir{irun(i)}];
end
Srev = [];
Srev.isubj = isubj;
Srev.irun = irun;
Srev.slist = dinfo;
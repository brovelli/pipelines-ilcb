function isval = cp_meg_review_gui(Sdb, opt)
% General GUI to review preprocessing parameters
% Last chance to modify channels / trials to remove
%-CREx180726

Srev = prepare_list(Sdb);
Srev.Sdb = Sdb;
Srev.opt = opt;
uiwait(msgbox({'\fontsize{12}Please confirm/change preprocessing parameters ';...
    'before source/connectivity analysis'}, 'Preprocessing review', 'help',...
    struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
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



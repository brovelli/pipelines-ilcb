function varargout = preproc_review(varargin)
% PREPROC_REVIEW MATLAB code for preproc_review.fig
%      PREPROC_REVIEW, by itself, creates a new PREPROC_REVIEW or raises the existing
%      singleton*.
%
%      H = PREPROC_REVIEW returns the handle to a new PREPROC_REVIEW or the handle to
%      the existing singleton*.
%
%      PREPROC_REVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROC_REVIEW.M with the given input arguments.
%
%      PREPROC_REVIEW('Property','Value',...) creates a new PREPROC_REVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preproc_review_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preproc_review_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preproc_review

% Last Modified by GUIDE v2.5 14-Aug-2018 02:10:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preproc_review_OpeningFcn, ...
                   'gui_OutputFcn',  @preproc_review_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before preproc_review is made visible.
function preproc_review_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preproc_review (see VARARGIN)

% Get the input param
Srev = varargin{1};
handles.list_subj.String = Srev.slist;
handles.Srev = Srev;
handles.change = 0;
handles.isel = 1;
handles.list_subj.Value = 1;
handles = disp_param(handles);
% Choose default command line output for preproc_review
handles.output = hObject;
set(handles.main, 'CloseRequestFcn', []);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preproc_review wait for user response (see UIRESUME)
uiwait(handles.main);


% --- Outputs from this function are returned to the command line.
function varargout = preproc_review_OutputFcn(hObject, eventdata, handles)  %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.change;  
delete(hObject)

% --- Executes on selection change in list_subj.
function list_subj_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to list_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_subj contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_subj

handles.isel = hObject.Value;
handles = disp_param(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function list_subj_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to list_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in but_ok.
function but_ok_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to but_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)
uiresume;

% --- Executes on button press in but_rs.
function but_rs_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to but_rs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Sdb, prevs] = set_new(handles, 'new_vis', 'sens');
handles.Srev.Sdb = cp_meg_rmsens_gui(Sdb, handles.Srev.opt.continuous);
handles = compare_aft(handles, 'sens', prevs);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in but_ra.
function but_ra_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to but_ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in but_rc.
function but_rc_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to but_rc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Sdb, prevc] = set_new(handles, 'new_ica', 'comp');
handles.Srev.Sdb = cp_meg_rmcomp_gui(Sdb);
handles = compare_aft(handles, 'comp', prevc);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in but_rt.
function but_rt_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to but_rt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Sdb, prevt] = set_new(handles, 'new_rmt', 'trials');
handles.Srev.Sdb = cp_meg_rmtrials_gui(Sdb);
handles = compare_aft(handles, 'trials', prevt);

% Update handles structure
guidata(hObject, handles);
%------------------------------------------------------------- Additionnal
% functions

function handles = compare_aft(handles, rmnam, prevrm)

isel = handles.isel;
i = handles.Srev.isubj(isel);
j = handles.Srev.irun(isel);
aftrm = handles.Srev.Sdb(i).meg.preproc.param_run{j}.rm.(rmnam);
if ~isstruct(prevrm)
    if ~isempty(setxor(aftrm, prevrm))
        handles.change = 1;
    end
else
   [cond, Nc] = get_names(prevrm);
    for k = 1 : Nc
        cnam = cond{k};
        if ~isempty(setxor(aftrm.(cnam), prevrm.(cnam)))
            handles.change = 1;
            break;
        end
    end 
end

function [Sdb, prevelm] = set_new(handles, newnam, rmnam)
Sdb = handles.Srev.Sdb;
isel = handles.isel;
i = handles.Srev.isubj(isel);
j = handles.Srev.irun(isel);
Sdb(i).meg.preproc.(newnam)(j) = 1;
prevelm = Sdb(i).meg.preproc.param_run{j}.rm.(rmnam);

function clist = parlist(clist)
if isempty(clist)
    clist = 'None';
    return;
end
if ~iscell(clist)
    if length(clist(1,:)) == 2
        clist = cellfun(@(x) add_cro(strjoint(strsplitt(x, ' '), '-')), clist, 'UniformOutput', 0); 
        clist = strjoint(clist, ' ; ');
    else
        % Assume numeric
        clist = cellstr(num2str(clist));
        clist = cellfun(@(x) strrep(x, ' ', ''), clist, 'UniformOutput', 0);
        clist = strjoint(clist, ', ');
    end
else
    clist = strjoint(clist, ' ; ');
end

function x = add_cro(x)
x = ['[', x, '] s'];
  
function clist = parlist_trials(Str)
[cond, Nc] = get_names(Str);
clist = cell(Nc, 1);
for k = 1 : Nc
    cnam = cond{k};
    if strcmp(cnam, 'allcond')
        cnamt = 'Same for all conditions: ';
    else
        cnamt = ['Condition ', cnam,': '];
    end
    clist{k} = [cnamt, parlist(Str.(cnam))];
end
clist = strjoint(clist, ' ; ');

function handles = disp_param(handles)
isel = handles.isel;
Srev = handles.Srev;
handles.panel.Title = Srev.slist{isel};
i = Srev.isubj(isel);
j = Srev.irun(isel);
Sdb = handles.Srev.Sdb;
Srm = Sdb(i).meg.preproc.param_run{j}.rm;
handles.lablist_rs.String = parlist(Srm.sens);
handles.lablist_ra.String = parlist(Srm.art);
handles.lablist_rc.String = parlist(Srm.comp);
handles.lablist_rt.String = parlist_trials(Srm.trials);

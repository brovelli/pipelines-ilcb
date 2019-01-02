function varargout = preproc_select(varargin)
% PREPROC_SELECT MATLAB code for preproc_select.fig
%      PREPROC_SELECT, by itself, creates a new PREPROC_SELECT or raises the existing
%      singleton*.
%
%      H = PREPROC_SELECT returns the handle to a new PREPROC_SELECT or the handle to
%      the existing singleton*.
%
%      PREPROC_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROC_SELECT.M with the given input arguments.
%
%      PREPROC_SELECT('Property','Value',...) creates a new PREPROC_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preproc_select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preproc_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preproc_select

% Last Modified by GUIDE v2.5 14-Dec-2018 17:02:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preproc_select_OpeningFcn, ...
                   'gui_OutputFcn',  @preproc_select_OutputFcn, ...
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


% --- Executes just before preproc_select is made visible.
function preproc_select_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preproc_select (see VARARGIN)

% Choose default command line output for preproc_select
handles.output = hObject;

% Get the input param
Sdisp = varargin{1};

Sdisp = prepare_listbox(Sdisp);

if isempty(Sdisp.dir)
    handles.push_opendir.Visible = 'off';
end
handles.dir = Sdisp.dir;
handles.list_good.String = Sdisp.good.sbox;
handles.list_good.Value = [];
handles.list_bad.String = Sdisp.bad.sbox;
handles.list_bad.Value = [];
handles.good.isel = [];
handles.good.clist = Sdisp.good.clist;
handles.bad.isel = [];
handles.bad.clist = Sdisp.bad.clist;
handles.title.String = Sdisp.title;
handles.title.BackgroundColor = [1 1 1];
handles.lab_good.BackgroundColor = [1 1 1];
handles.lab_bad.BackgroundColor = [1 1 1];
handles.lab_click.BackgroundColor = [1 1 1];
handles.badrun = 0;
set(handles.main, 'CloseRequestFcn', []);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preproc_select wait for user response (see UIRESUME)
uiwait(handles.main);

% --- Outputs from this function are returned to the command line.
function varargout = preproc_select_OutputFcn(hObject, eventdata, handles)  %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.bad.clist;
varargout{2} = handles.badrun;
delete(hObject)

% --- Executes on button press in push_bad.
function push_bad_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to push_bad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isel = handles.good.isel;
if isempty(isel)
    return
end
csel = handles.good.clist(isel);
ckeep = setxor(handles.good.clist, csel);
handles.good.clist = ckeep;
handles.bad.clist = [handles.bad.clist ; csel];
handles = list_update(handles);
guidata(hObject,handles)

% --- Executes on button press in push_good.
function push_good_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to push_good (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isel = handles.bad.isel;
if isempty(isel)
    return
end
csel = handles.bad.clist(isel);
ckeep = setxor(handles.bad.clist, csel);
handles.bad.clist = ckeep;
handles.good.clist = [handles.good.clist ; csel];
handles = list_update(handles);
guidata(hObject,handles)

% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to push_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)
uiresume;

% --- Executes on selection change in list_good.
function list_good_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to list_good (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_good contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_good
handles.good.isel = hObject.Value;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function list_good_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to list_good (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_bad.
function list_bad_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to list_bad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_bad contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_bad
handles.bad.isel = hObject.Value;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function list_bad_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to list_bad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in push_opendir.
function push_opendir_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to push_opendir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pfig = handles.dir;
iso = 1;
if ispc
    winopen(pfig)
elseif ismac
    try
        system(['open ', pfig, ' &']);
    catch
        iso = 0;
    end    
end
if ~iso
    smsg = {'Impossible to open the processing directory with Matlab';
        'Figures can be found here: ';
        pfig};
    uiwait(msg_box(smsg));
end

% --- Executes on button press in push_badrun.
function push_badrun_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to push_badrun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.badrun = 1;
guidata(hObject,handles)
uiresume;
%------------------------ Additionnal functions

function Sdisp = prepare_listbox(Sdisp)
Sdisp.good = strlist(Sdisp.good, [0 0.7 0]);
Sdisp.bad = strlist(Sdisp.bad, [0.75 0 0]);

function S = strlist(S, col)

clist = S.clist;
if isempty(clist)
    S.sbox = '';
    return
end

col = rgb2hex(col);
pre = ['<html><font color="FFFFFF">x </font><font color="', col,'">'];
post = '</font></p></html>';

if ~iscell(clist) && isnumeric(clist(1))
    slist = sort(clist);
    clist = arrayfun(@num2str, slist, 'UniformOutput', 0);
else
    if strcmp(clist{1}(1), 'A')
        numc = cellfun(@(x) str2double(x(2:end)), clist, 'UniformOutput', 1);
        [~, inds] = sort(numc);
        clist = clist(inds);
    end
    slist = clist;
end

ndm = cellfun(@length, clist, 'UniformOutput', 1);
N = max(ndm);

S.sbox = cellfun(@(x) [pre, blanks(N-length(x)), x, post], clist, 'UniformOutput', 0);
S.clist = slist;
% See 
% https://fr.mathworks.com/matlabcentral/answers/153064-change-color-of-each-individual-string-in-a-listbox
function hex = rgb2hex(col)
col = round(col*255);
hex = reshape(dec2hex(col(1, :), 2)', 1, 6);

function handles = list_update(handles)

handles = prepare_listbox(handles);
handles.good.isel = [];
handles.bad.isel = [];
handles.list_good.String = handles.good.sbox;
handles.list_good.Value = [];
handles.list_bad.String = handles.bad.sbox;
handles.list_bad.Value = [];

pause(0.05)

% --- Executes when user attempts to close main_param.
% function main_CloseRequestFcn(hObject, eventdata, handles)  %#ok
% % hObject    handle to main_param (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % varargout{1} = handles.list_bad.UserData;  %#ok
% % delete(hObject)
% set(handles.main, 'visible', 'off')
% Hint: delete(hObject) closes the figure

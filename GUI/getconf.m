function answer = getconf(varargin)

if nargin == 0  % LAUNCH GUI

  fig = openfig(mfilename,'reuse');

  % Use system color scheme for figure:
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
  
  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  guidata(fig, handles);
  initial_dir = pwd;
  % Populate the listbox
  load_listbox(initial_dir,handles,0)
  
  uiwait(fig);
  handles=guidata(fig);
  %answer = listbox1_Callback(handles.listbox1, [], handles, varargin)
  % might have returned because the window was deleted using
  % the close box - in that case, return 'cancel' as the answer, and
  % don't bother deleting the window!
  if ~ishandle(fig)
      %      handles
      answer = 'cancel';
      
  else
      % otherwise, we got here because the user pushed one of the two buttons.
      % retrieve the latest copy of the 'handles' struct, and return the answer.
      % Also, we need to delete the window.
      answer = handles.answer;
      delete(fig);
  end
  
  if nargout > 0
      varargout{1} = fig;
  end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
  try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
  catch
    disp(lasterr);
  end

else  % LAUNCH GUI

  fig = openfig(mfilename,'reuse');

  % Use system color scheme for figure:
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
  
  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  set(handles.popupmenu1,'Value',varargin{1})
  guidata(fig, handles);
  initial_dir = pwd;
  % Populate the listbox
  load_listbox(initial_dir,handles,varargin{1})
  
  uiwait(fig);
  handles=guidata(fig);
  %answer = listbox1_Callback(handles.listbox1, [], handles, varargin)
  % might have returned because the window was deleted using
  % the close box - in that case, return 'cancel' as the answer, and
  % don't bother deleting the window!
  if ~ishandle(fig)
      %      handles
      answer = 'cancel';
      
  else
      % otherwise, we got here because the user pushed one of the two buttons.
      % retrieve the latest copy of the 'handles' struct, and return the answer.
      % Also, we need to delete the window.
      answer = handles.answer;
      delete(fig);
  end
end




  
% ------------------------------------------------------------
% Callback for list box - open .fig with guide, otherwise use open
% ------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)

if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.listbox1,'Value');
    file_list = get(handles.listbox1,'String');	
    filename = file_list{index_selected};
    if  handles.is_dir(handles.sorted_index(index_selected))
        cd (filename)
        load_listbox(pwd,handles,0)
    else
        handles.answer=filename;
        guidata(h, handles);
        uiresume(handles.figure1);
    end
else
    list_entries = get(handles.listbox1,'String');
    index_selected = get(handles.listbox1,'Value');
    set(handles.edit1,'String',list_entries{index_selected(1)});
    guidata(h, handles);
end

% ------------------------------------------------------------
% Read the current directory and sort the names
% ------------------------------------------------------------
function load_listbox(dir_path,handles,v)
cd (dir_path)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
[n m]=size(sorted_names);
conf_names=sorted_names(1);
conf_index=sorted_index(1);
handles.is_dir = [dir_struct.isdir];
guidata(handles.figure1,handles);

if (v>0 & v<9)
    j=1;
    switch v
    case 1
        file_type = '.conf';
    case 2
        file_type = '.tie';
    case 3
        file_type = '.gap';
    case 4
        file_type = '.txt';
    case 5
        file_type = '.new';
    case 6
        file_type = '.eqn';
    case 7
        file_type = '.offset';
    case 8
        file_type = '.mcdtie';
    end
    
    for i=1:n
        if  handles.is_dir(sorted_index(i))
            conf_names(j)=sorted_names(i);
            conf_index(j)=sorted_index(i);
            j=j+1;
        else
            [path,name,ext] = fileparts(sorted_names{i});  
            switch  ext
            case file_type
                conf_names(j)=sorted_names(i);
                conf_index(j)=sorted_index(i);
                j=j+1;
            end
        end
    end
else
    conf_names=sorted_names;
    conf_index=sorted_index;
end

handles.file_names = conf_names;

handles.sorted_index = [conf_index];

set(handles.listbox1,'String',handles.file_names,'Value',1)
set(handles.text1,'String',pwd)
guidata(handles.figure1,handles)


% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.


handles.answer=get(handles.edit1,'String');
guidata(h, handles);
uiresume(handles.figure1);


% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit1.

% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.

handles.answer='cancel';
guidata(h, handles);
uiresume(handles.figure1);

% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu1.

load_listbox(pwd,handles,get(handles.popupmenu1,'Value'));

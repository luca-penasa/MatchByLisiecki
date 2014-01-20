function answer = getspeed(varargin)

if ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

  try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
  catch
    disp(lasterr);
  end

else

  fig = openfig(mfilename,'reuse');
 	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

 % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  handles.num=varargin{1};  
  handles.den=varargin{2};
  guidata(fig, handles);
  s=speeds2str(handles);
  set(handles.edit1,'String',s);
  handles.speed=s;
  handles.s=s;
  guidata(fig, handles);
  
  uiwait(fig);
  handles=guidata(fig);
  
  if ~ishandle(fig)
      answer = s;
  else
  	  % otherwise, we got here because the user pushed one of the two buttons.
	  % retrieve the latest copy of the 'handles' struct, and return the answer.
	  % Also, we need to delete the window.
      answer = handles.speed;
	  delete(fig);
  end

   
end



% --------------------------------------------------------------------
function answer = speeds2str(handles)
    
[m n]=size(handles.num);
s='';
for i=1:n
    s(end+1)=num2str(handles.num(i));
    s(end+1)=':';
    s(end+1)=num2str(handles.den(i));
    if(i~=n)
        s(end+1)=',';
    end
end
answer=s;

% --------------------------------------------------------------------
function answer = verify(h, eventdata, handles)

s='';
s=get(handles.edit1,'String');
if(strcmp(s,''))
    answer=0;
else
    ind=find(s==':');
    [m n]=size(ind);
    if (n>0 & max(size(s))>=(4*n-2) & min(ind)>1 & max(ind)<max(size(s)))
        answer=1;
    else
        answer=0;
    end
end
        
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit1.
good=verify(h, eventdata, handles);
if good
    handles.speed=get(handles.edit1,'String');
    guidata(h, handles);
else
    set(handles.edit1,'String','Must contain >1 ratio (e.g. 1:2,1:1)');
end

% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
uiresume(handles.figure1);



% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.
set(handles.edit1,'String',handles.s);
guidata(h, handles);

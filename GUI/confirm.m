function answer = confirm(varargin)
%function varargout = confirm(varargin)
% CONFIRM Application M-file for confirm.fig
%    FIG = CONFIRM launch confirm GUI.
%    CONFIRM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 24-Nov-2003 15:35:25

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);
    handles=guidata(fig);
    % might have returned because the window was deleted using
    % the close box - in that case, return 'cancel' as the answer, and
    % don't bother deleting the window!
    if ~ishandle(fig)
        %      handles
        answer = 0;
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

else
    fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    v=varargin{1};
    if(v==1)
        set(handles.text1,'String','Exit Match configuration interface?');
    elseif(v==2)
        set(handles.text1,'String','Delete all tie points?');
    elseif(v==3)
        set(handles.text1,'String','Delete all user-defined gaps?');
    elseif(v==4)
        set(handles.text1,'String','Load gaps for signal or target?');
        set(handles.cancel_button,'String','Signal');   
        set(handles.yes_button,'String','Target');   
    elseif(v==5)
        set(handles.text1,'String','Apply MCD conversion and exit?');
    elseif(v==6)
        set(handles.text1,'String','Exit MCD conversion but do not apply?');
    elseif(v==7)
        set(handles.text1,'String','MCD tie points for linear interpolation:');
        set(handles.yes_button,'String','Hole 1 + User-Defined','FontSize',7);   
        set(handles.cancel_button,'String','User-Defined Only','FontSize',7);   
    end
	guidata(fig, handles);
    
    uiwait(fig);
    
    % might have returned because the window was deleted using
    % the close box - in that case, return 'cancel' as the answer, and
    % don't bother deleting the window!
    if ~ishandle(fig)
        %      handles
        answer = 0;
        
    else
        % otherwise, we got here because the user pushed one of the two buttons.
        % retrieve the latest copy of the 'handles' struct, and return the answer.
        % Also, we need to delete the window.
        handles=guidata(fig);
        answer = handles.answer;
        delete(fig);
    end
    
    
    if nargout > 0
        varargout{1} = fig;
    end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = yes_button_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.yes_button.
handles.answer=1;
guidata(h, handles);
uiresume(handles.figure1);


% --------------------------------------------------------------------
function varargout = cancel_button_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel_button.
handles.answer=0;
guidata(h, handles);
uiresume(handles.figure1);


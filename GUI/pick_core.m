function answer = pick_core(varargin)
% PICK_CORE Application M-file for pick_core.fig
%    FIG = PICK_CORE launch pick_core GUI.
%    PICK_CORE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 28-Jul-2003 13:02:12

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    
    uiwait(fig);
    handles=guidata(fig);
    %answer = listbox1_Callback(handles.listbox1, [], handles, varargin)
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
    v=varargin{1}
    c2=['core ' num2str(v(1)) ' which extends from ' num2str(v(2)) ' to ' num2str(v(3))];
    c1=['core ' num2str(v(4)) ' which extends from ' num2str(v(5)) ' to ' num2str(v(6))];
    question=['Place in ' c1  ' or in ' c2 '?']; 
    set(handles.text1,'String',question);
    c1=['Core ' num2str(v(4))];
    set(handles.pushbutton1,'String',c1);
    c2=['Core ' num2str(v(1))];
    set(handles.pushbutton2,'String',c2);
	guidata(fig, handles);
    
    uiwait(fig);
    handles=guidata(fig);
    %answer = listbox1_Callback(handles.listbox1, [], handles, varargin)
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
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

handles.answer=1;
guidata(h, handles);
uiresume(handles.figure1);


% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.

handles.answer=2;
guidata(h, handles);
uiresume(handles.figure1);

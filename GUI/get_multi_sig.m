function files = get_multi_sig(varargin)
% GET_MULTI_SIG Application M-file for get_multi_sig.fig
%    FIG = GET_MULTI_SIG launch get_multi_sig GUI.
%    GET_MULTI_SIG('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 01-Aug-2003 11:44:58

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    uiwait(fig);
    handles=guidata(fig);
    % might have returned because the window was deleted using
    % the close box - in that case, return 'cancel' as the answer, and
    % don't bother deleting the window!
    if ~ishandle(fig)
        %      handles
        files = [];
        
    else
        % otherwise, we got here because the user pushed one of the two buttons.
        % retrieve the latest copy of the 'handles' struct, and return the answer.
        % Also, we need to delete the window.
        answer = handles.answer;
        if(answer==1)
            s2=get(handles.sig2,'String');
            t2=get(handles.targ2,'String');
            s3=get(handles.sig3,'String');
            t3=get(handles.targ3,'String');
            s4=get(handles.sig4,'String');
            t4=get(handles.targ4,'String');
            if((~isempty(s2) & isempty(t2)) | (isempty(s2) & ~isempty(t2)))
                files=[];
                disp('Error: Must enter both signal2 and target2')
            elseif((~isempty(s3) & isempty(t3)) | (isempty(s3) & ~isempty(t3)))
                files=[];
                disp('Error: Must enter both signal 3 and target 3')
            elseif((~isempty(s4) & isempty(t4)) | (isempty(s4) & ~isempty(t4)))
                files=[];
                disp('Error: Must enter both signal 4 and target4')
            else
                files=char(s2, t2, s3, t3, s4, t4);
            end
        else
            files = [];
        end
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
function varargout = sig2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig2.

% --------------------------------------------------------------------
function varargout = targ2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.targ2.

% --------------------------------------------------------------------
function varargout = sig3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig3.

% --------------------------------------------------------------------
function varargout = targ3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.targ3.

% --------------------------------------------------------------------
function varargout = sig4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig4.

% --------------------------------------------------------------------
function varargout = targ4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.targ4.


% --------------------------------------------------------------------
function varargout = accept_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.accept.
handles.answer=1;
guidata(h, handles);
uiresume(handles.figure1);


% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel.
handles.answer=0;
guidata(h, handles);
uiresume(handles.figure1);

% --------------------------------------------------------------------
function varargout = browse_s2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.sig2,'String',sig);
end
    
% --------------------------------------------------------------------
function varargout = browse_t2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.targ2,'String',sig);
end
    
% --------------------------------------------------------------------
function varargout = browse_s3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.sig3,'String',sig);
end
    
% --------------------------------------------------------------------
function varargout = browse_t3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.targ3,'String',sig);
end
    
% --------------------------------------------------------------------
function varargout = browse_s4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.sig4,'String',sig);
end
    
% --------------------------------------------------------------------
function varargout = browse_t4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.targ4,'String',sig);
end
    

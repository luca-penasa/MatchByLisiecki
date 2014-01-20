function answer = insert_mcdtie(varargin)
% INSERT_MCDTIE Application M-file for insert_mcdtie.fig
%    FIG = INSERT_MCDTIE launch insert_mcdtie GUI.
%    INSERT_MCDTIE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 06-Jun-2006 18:46:29

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    guidata(fig, handles);

    set(handles.radio_ccd,'String','MCD to CMCD (enter any 2 values)');
    set(handles.radio_mbsf,'String','MBSF to CMCD');
    set(handles.text5,'String','CMCD');
    set(handles.text6,'String','MCD');
    set(handles.text7,'String','CMCD');

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
function varargout = radio_mbsf_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_mbsf.
if(get(handles.radio_mbsf,'Value')==1)
    set(handles.radio_ccd,'Value',0);
else
    set(handles.radio_ccd,'Value',1);
end


% --------------------------------------------------------------------
function varargout = radio_ccd_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccd.
if(get(handles.radio_ccd,'Value')==1)
    set(handles.radio_mbsf,'Value',0);
else
    set(handles.radio_mbsf,'Value',1);
end


% --------------------------------------------------------------------
function varargout = edit_hole_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_hole.
if(get(handles.radio_mbsf,'Value')==0)
    set(handles.radio_mbsf,'Value',1);
    set(handles.radio_ccd,'Value',0);
end


% --------------------------------------------------------------------
function varargout = edit_mbsf_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_mbsf.
if(get(handles.radio_mbsf,'Value')==0)
    set(handles.radio_mbsf,'Value',1);
    set(handles.radio_ccd,'Value',0);
end


% --------------------------------------------------------------------
function varargout = edit_mcd1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_mcd1.
if(get(handles.radio_mbsf,'Value')==0)
    set(handles.radio_mbsf,'Value',1);
    set(handles.radio_ccd,'Value',0);
end


% --------------------------------------------------------------------
function varargout = edit_ccd_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_ccd.
if(get(handles.radio_ccd,'Value')==0)
    set(handles.radio_mbsf,'Value',0);
    set(handles.radio_ccd,'Value',1);
end


% --------------------------------------------------------------------
function varargout = edit_mcd2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_mcd2.
if(get(handles.radio_ccd,'Value')==0)
    set(handles.radio_mbsf,'Value',0);
    set(handles.radio_ccd,'Value',1);
end


% --------------------------------------------------------------------
function varargout = edit_offset_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_offset.
if(get(handles.radio_ccd,'Value')==0)
    set(handles.radio_mbsf,'Value',0);
    set(handles.radio_ccd,'Value',1);
end


% --------------------------------------------------------------------
function varargout = push_ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_ok.
if(get(handles.radio_mbsf,'Value')==1)
    f1=get(handles.edit_hole,'String');
    f2=get(handles.edit_mbsf,'String');
    f3=get(handles.edit_mcd1,'String');
    try
        if(~isempty(f1) & (isempty(str2num(f1)) | f1=='i' | f1=='j'))
            f1=num2str(double(f1));
            typept=0;
        else
            typept=1;
        end
        if(~isempty(str2num(f1)) & ~isempty(str2num(f2)) & ~isempty(str2num(f3)))
            handles.answer=[str2num(f1) str2num(f2) str2num(f3) typept];
            guidata(h, handles);
            %delete(handles.mcdfigure);
            uiresume(handles.insertfig);
        else
            disp('Error 1: Must enter all 3 values for hole/mbsf/mcd')
        end
    catch
        disp('Error 2: Must enter all 3 values for hole/mbsf/mcd')
    end
    
else
    f1=get(handles.edit_ccd,'String');
    f2=get(handles.edit_mcd2,'String');
    f3=get(handles.edit_offset,'String');
    b=[0 0 0];
    % a(1)=CCD a(2)=MCD a(3)=offset=a(1)-a(2)
    a=[NaN NaN NaN 2];
    if(~isempty(f1))
        b(1)=~isempty(str2num(f1));
    end
    if(~isempty(f2))
        b(2)=~isempty(str2num(f2));
    end
    if(~isempty(f3))
        b(3)=~isempty(str2num(f3));
    end
    
    if(sum(b)>=2)
        for i=1:3
            if(b(i)==1)
                a(i)=str2num(eval(['f' num2str(i)]));
            end
        end        
        if(isnan(a(1)))
            a(1)=a(2)+a(3);
        end
        if(isnan(a(2)))
            a(2)=a(1)-a(3);
        end
        if(isnan(a(3)))
            a(3)=a(1)-a(2);
        end
        
        handles.answer=a;
        guidata(h, handles);
        %delete(handles.mcdfigure);
        uiresume(handles.insertfig);
    else
        disp('Error: Must fill in 2 of the 3 boxes for ccd/mcd/offset')
    end
end

% --------------------------------------------------------------------
function varargout = push_cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_cancel.
handles.answer=[NaN NaN NaN NaN];
guidata(h, handles);
%delete(handles.mcdfigure);
uiresume(handles.insertfig);

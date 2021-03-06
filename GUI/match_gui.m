% This program generates configuration files for and displays the results
% of Match 2.0,a signal correlation program designed for aligning paleoclimate
% records. A user's manual is available. Questions and update requests can
% be directed to Lorraine Lisiecki, zogalum@alum.mit.edu.

function varargout = match_gui(varargin)
% Match_GUI Application M-file for match_gui.fig
%    FIG = Match_GUI launch match_gui GUI.
%    Match_GUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 05-Aug-2003 15:08:17

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    disp(' ')
    disp('Welcome to the Match 2.0 graphical user interface.')
    disp('Created by Lorraine Lisiecki, 2003, with funding from')
    disp('a Schlanger Ocean Drilling Fellowship (NSF-USSSP).')
    disp(' ')

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.nties=1;
    handles.ties=[];
    handles.ngaps=1;
    handles.gaps=[];
    handles.stop=0;
    handles.xc1=0;
    handles.yc1=0;
    handles.xc2=0;
    handles.yc2=0;
    handles.speed_num=[1,2,1,3,2,3,4,1,5,4,3,5,2,5,3];
    handles.speed_den=[3,5,2,5,3,4,5,1,4,3,2,3,1,2,1];
	handles.target_speed='1:1';
    handles.oldsig=get(handles.edit_sig,'String');
    handles.oldtarg=get(handles.edit_targ,'String');
    handles.moresig=[];
    handles.sig_fig=0;
    handles.match_fig=0;
    handles.results_axes=[];
    guidata(fig, handles);
    dummy=speeds2str(handles);

    
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
function varargout = edit_sig_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_sig.
% --------------------------------------------------------------------
% Reads in signal from text box (handles.edit_sig) and plots file contents;
% also sets default signal start and end points.

handles=guidata(handles.figure1);

f1=get(handles.edit_sig,'String');
n=0;
try
    s1=load(f1);
    [m n]=size(s1);
end

% Test that file has the required 2 or 3 columns
if(n==3 | n==2)
    xc=n-1;
    yc=n;
    handles.xc1=xc;
    handles.yc1=yc;
    st=min(s1(:,xc));
    e=max(s1(:,xc));
    set(handles.slider1,'min',st,'max',e,'Value',st,'sliderstep',[0 0]);
    handles.s1max=max(s1(:,yc));
    handles.x1max=e;
    guidata(gcbo,handles);
    axes(handles.axes1)
    % delete any pre-existing tie points
    if(handles.nties>1)
        push_del_all_Callback(h, eventdata, handles,1)
    end
    handles=guidata(handles.figure1);
    
    % delete any pre-existing signal gaps
    if(~isempty(handles.gaps))
        ind=find(handles.gaps(:,2)==1);
        if(~isempty(ind))
            handles.gaps(ind,:)=[];
            guidata(gcbo,handles);
        end
    end
    
    % find gaps in signal denoted by NaN
    [m,p]=find(isnan(s1));
    j=size(m);
    % find gaps in signal denoted by changes in core number
    if(yc==3 & isempty(m))
        [m,p]=find(s1(1:end-1,1)~=s1(2:end,1));
        j=size(m);
        for i=1:j
            s1(m(i)+i:end+1,:)=[NaN NaN NaN;s1(m(i)+i:end,:)];
            m(i)=m(i)+i-1;
        end
    elseif(~isempty(m))
        m=m-1
    end
    
    % plot signal
    hold off;
    plot(s1(:,xc),s1(:,yc),'m');
    hold on;
    % plot red dots at end of each core
    for i=1:j
        plot(s1(m(i),xc),s1(m(i),yc),'.','Color',[1 .01 .01],'MarkerSize',9);
    end
    xlim(handles.axes1,[st e]);
    flipx_Callback(h, eventdata, handles)
    flipy_Callback(h, eventdata, handles)
    b=axis;
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim1,'String',c);
    set(gca,'DefaulttextButtonDownFcn','match_gui(''text1_ButtondownFcn'',gcbo,[],guidata(gcbo))');
    guidata(gcbo,handles);
    
    if(isempty(varargin))
        % set begin and end points to start and end of signal
        set(handles.edit_beg1,'String',num2str(st));
        set(handles.edit_end1,'String',num2str(e));
        delete(findobj('Marker','>','Parent',handles.axes1));
        plot(st,s1(1,yc),'g>','MarkerSize',7,'MarkerFaceColor','g');
        delete(findobj('Marker','<','Parent',handles.axes1));
        plot(e,s1(end,yc),'g<','MarkerSize',7,'MarkerFaceColor','g');
    else
        % plot pre-existing begin and end points if signal has not changed
        e=str2num(get(handles.edit_end1,'String'));
        x=find(s1(:,xc)<=e);
        st=str2num(get(handles.edit_beg1,'String'));
        plot(e,s1(x(end),yc),'g<','MarkerSize',7,'MarkerFaceColor','g');
        x=find(s1(:,xc)>=st);
        plot(st,s1(x(1),yc),'g>','MarkerSize',7,'MarkerFaceColor','g');
    end
    
    handles.oldsig=get(handles.edit_sig,'String');
    guidata(gcbo,handles);
    
else
    % If file does not have 2 or 3 columns revert to previous signal
    disp('Error: unrecognized series1 format')
    set(handles.edit_sig,'String',handles.oldsig);
end

% --------------------------------------------------------------------
function varargout = edit_targ_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_targ.
% --------------------------------------------------------------------
% Reads in target from text box (handles.edit_targ) and plots file contents;
% also sets default signal start and end points.

handles=guidata(handles.figure1);
ntie=handles.nties;
f2=get(handles.edit_targ,'String');
n=0;
try
    s2=load(f2);
    [m n]=size(s2);
end

if(n==3 | n==2)
    xc=n-1;
    yc=n;
    handles.xc2=xc;
    handles.yc2=yc;
    st=min(s2(:,xc));
    e=max(s2(:,xc));
    set(handles.slider2,'min',st,'max',e,'Value',st,'sliderstep',[0 0]);
    handles.s2max=max(s2(:,yc));
    handles.x2max=e;
    guidata(gcbo,handles);
    axes(handles.axes2)
    %delete pre-existing tie points
    if(handles.nties>1)
        push_del_all_Callback(h, eventdata, handles,1)
    end
    handles=guidata(handles.figure1);
    %delete pre-existing gaps in target
    if(~isempty(handles.gaps))
        ind=find(handles.gaps(:,2)==2);
        if(~isempty(ind))
            handles.gaps(ind,:)=[];
            guidata(gcbo,handles);
       end
    end
    % find gaps in signal denoted by NaN 
    [m,p]=find(isnan(s2));
    j=size(m);
    % find gaps in signal denoted by changes in core number
    if(yc==3 & isempty(m))
        [m,p]=find(s2(1:end-1,1)~=s2(2:end,1));
        j=size(m);
        for i=1:j
            s2(m(i)+i:end+1,:)=[NaN NaN NaN;s2(m(i)+i:end,:)];
            m(i)=m(i)+i-1;
        end
    elseif(~isempty(m))
        m=m-1;
    end
    
    %plot target
    hold off;
    plot(s2(:,xc),s2(:,yc),'Color',[.35 .35 1]);
    hold on;
    %plot red dots at end of each core
    for i=1:j
        plot(s2(m(i),xc),s2(m(i),yc),'.','Color',[1 .01 .01],'MarkerSize',9);
    end
    xlim(handles.axes2,[st e]);
    flipx2_Callback(h, eventdata, handles)
    flipy2_Callback(h, eventdata, handles)
    b=axis;
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim2,'String',c);
    set(gca,'DefaulttextButtonDownFcn','match_gui(''text2_ButtondownFcn'',gcbo,[],guidata(gcbo))');
    
    if(isempty(varargin))
        % set begin and end points to start and end of signal if new target
        set(handles.edit_beg2,'String',num2str(st));
        set(handles.edit_end2,'String',num2str(e));
        delete(findobj('Marker','>','Parent',handles.axes2));
        plot(st,s2(1,yc),'g>','MarkerSize',7,'MarkerFaceColor','g');
        delete(findobj('Marker','<','Parent',handles.axes2));
        plot(e,s2(end,yc),'g<','MarkerSize',7,'MarkerFaceColor','g');
    else
        % plot pre-existing begin and end points if target has not changed
        e=str2num(get(handles.edit_end2,'String'));
        x=find(s2(:,xc)<=e);
        st=str2num(get(handles.edit_beg2,'String'));
        plot(e,s2(x(end),yc),'g<','MarkerSize',7,'MarkerFaceColor','g');
        x=find(s2(:,xc)>=st);
        plot(st,s2(x(1),yc),'g>','MarkerSize',7,'MarkerFaceColor','g');
    end
    handles.oldtarg=get(handles.edit_targ,'String');
    guidata(gcbo,handles);
    
else
    % If file does not have 2 or 3 columns revert to previous signal
    disp('Error: unrecognized series2 format')
    set(handles.edit_targ,'String',handles.oldtarg);
end

% --------------------------------------------------------------------
function varargout = load_button_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.load_button.
% --------------------------------------------------------------------
% Reads in pre-existing configuration file

% Get configuration filename
fconf=getconf(1);

if(~strcmp(fconf,'cancel'))
    tiefile=[];
    gapfile1=[];
    gapfile2=[];
    handles.gaps=[];
    handles.ngaps=1;
    % Read file
    a=textread(fconf,'%s');
    [m n]=size(a);
    % Store parameter values from file
    for i=1:m-1
        switch a{i}
        case 'series1',set(handles.edit_sig,'String',a{i+1});
            j=1;
            f=[];
            while(~strcmp('begin1',a{i+j+1}))
                j=j+1;
            end
            if(j>1)
                set(handles.multi_sig,'Value',1);
                f{1:2:2*(j-1)}=a{i+2:i+j};
            end
        case 'begin1',set(handles.edit_beg1,'String',a{i+1});
        case 'end1', set(handles.edit_end1,'String',a{i+1});
        case 'numintervals1', set(handles.editint,'String',a{i+1});
        case 'series2',set(handles.edit_targ,'String',a{i+1});
            if(j>1)
                f{2:2:2*(j-1)}=a{i+2:i+j};
                handles.moresig=char(f);
            end
        case 'begin2',set(handles.edit_beg2,'String',a{i+1});
        case 'end2',set(handles.edit_end2,'String',a{i+1});
        case 'numintervals2', set(handles.edit_nint2,'String',a{i+1});
        case 'nomatch', set(handles.end_pen,'String',a{i+1});
        case 'speedchange', set(handles.spch_pen,'String',a{i+1});
        case 'speedpenalty', set(handles.speed_pen,'String',a{i+1});
        case 'tiepenalty',  set(handles.tie_pen,'String',a{i+1});
        case 'gappenalty', set(handles.gap_pen,'String',a{i+1});
        case 'speeds', speeds=a{i+1};
        case 'tiefile', tiefile=a{i+1};
        case 'series1gaps', gapfile1=a{i+1};
        case 'series2gaps', gapfile2=a{i+1};
        case 'targetspeed', handles.target_speed=a{i+1};
        end
    end
    set(handles.edit_out,'String',fconf(1:end-5));
    guidata(gcbo,handles);
    
    % Load and plot signal and target
    edit_sig_Callback(h, eventdata, handles, 1);
    edit_targ_Callback(h, eventdata, handles, 1);
    guidata(gcbo,handles);

    % Load and plot tie points
    if(~isempty(tiefile))
        t=load(tiefile);
        [m n]=size(t);
        ind=[1:m]';
        if(n==2)
            handles.ties=[ind 0*ind t(:,1) 0*ind t(:,2)];
        else
            handles.ties=[ind t];
        end
        handles.nties=m+1;
        guidata(gcbo,handles);
        plot_ties(h, eventdata, handles,varargin);
    else
        handles.ties=[];
        handles.nties=1;
        guidata(gcbo,handles);
    end

    % Load and plot gaps from gap files
    for ax=1:2
        g=[];
        if(ax==1 & ~isempty(gapfile1))
            f1=get(handles.edit_sig,'String');
            g=load(gapfile1);
        end
        if(ax==2 & ~isempty(gapfile2))
            f1=get(handles.edit_targ,'String');
            g=load(gapfile2);
        end
        if(~isempty(g))
            try
                s1=load(f1);
                [m n1]=size(s1);
            end
            [m n]=size(g);
            
            for i=1:m
                x=g(i,2:3);        
                if(n1==2)
                    if(isnan(x(1)))
                        x(1)=min(s1(:,1));
                    end
                    if(isnan(x(2)))
                        x(2)=max(s1(:,1));
                    end
                else
                    if(isnan(x(1)))
                        x(1)=min(s1(find(s1(:,1)==g(i,1)),n1-1));
                    end
                    if(isnan(x(2)))
                        x(2)=max(s1(find(s1(:,1)==g(i,1)),n1-1));
                    end
                end
                handles.gaps(end+1,:)=[handles.ngaps ax g(i,1) x(1) x(2)];
                handles.ngaps=handles.ngaps+1;
            end
        end
        guidata(gcbo,handles);
        plot_ties(h, eventdata, handles, 1);
    end
    guidata(gcbo,handles);
 
    % Uncomment the following section to automatically plot multiple signals    
    %     if(get(handles.multi_sig,'Value')==1)
    %         plot_moresig(h, eventdata, handles);
    %     end


    %set match speeds
    ind=find(speeds(:)==':');
    ind2=find(speeds(:)==',' | speeds(:)==';' | speeds(:)=='[' | speeds(:)=='(');
    speeds(ind2)=' ';
    [m n]=size(ind);
    num=[];
    den=[];
    num(1)=sscanf(speeds(1:ind(1)-1),'%d');
    for i=1:m-1
        b=sscanf(speeds(ind(i)+1:ind(i+1)-1),'%d %d');
        den(end+1)=b(1);
        num(end+1)=b(2);
    end
    den(end+1)=sscanf(speeds(ind(m)+1:end),'%d');
    
    handles.speed_num=num;
    handles.speed_den=den;
    guidata(gcbo,handles);
    
    speeds=speeds2str(handles);
end


% --------------------------------------------------------------------
function varargout = editint_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.editint.
% Text box for number of matching intervals to use for signal

% --------------------------------------------------------------------
function varargout = edit_nint2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_nint2.
% Text box for number of matching intervals to use for target

% --------------------------------------------------------------------
function varargout = edit_out_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_out.
% --------------------------------------------------------------------
% Contains filename used as name base for configuration, tie, gap, match, 
% and log files

% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
% --------------------------------------------------------------------
% Generate configuration file

% Get values stored in text boxes and handles
str1=get(handles.edit_sig,'String');
str2=get(handles.edit_targ,'String');
nint1=str2num(get(handles.editint,'String'));
nint2=str2num(get(handles.edit_nint2,'String'));
nmpen=str2num(get(handles.end_pen,'String'));
sppen=str2num(get(handles.speed_pen,'String'));
spchpen=str2num(get(handles.spch_pen,'String'));
tiepen=str2num(get(handles.tie_pen,'String'));
gappen=str2num(get(handles.gap_pen,'String'));

labels=char('nomatch','speedchange','speedpenalty','gappenalty',...
    'tiepenalty','begin1','end1','numintervals1','begin2','end2',...
    'numintervals2','series1','series2','tiefile','matchfile','logfile');

% Read in signal and target
f1=str1;
f2=str2;
f1=get(handles.edit_sig,'String');
n1=0;
try
    s1=load(f1);
    [m n1]=size(s1);
    xc1=n1-1;
    yc1=n1;
end

n2=0;
try
    s2=load(f2);
    [m n2]=size(s2);
    xc2=n2-1;
    yc2=n2;
end

if(~(n1==2 | n1==3))
    disp('Error: Invalid signal- configuration file not written')
elseif(~(n2==2 | n2==3))
    disp('Error: Invalid target- configuration file not written')
else
    % Set begin and end to default if text boxes are empty
    if(isempty(get(handles.edit_end1,'String')))
        e1=s1(end,xc1);
    else
        e1=str2num(get(handles.edit_end1,'String'));
    end
    if(isempty(get(handles.edit_end2,'String')))
        e2=s2(end,xc2);
    else
        e2=str2num(get(handles.edit_end2,'String'));
    end
    
    if(isempty(get(handles.edit_beg1,'String')))
        st1=s1(1,xc1);
    else
        st1=str2num(get(handles.edit_beg1,'String'));
    end
    if(isempty(get(handles.edit_beg2,'String')))
        st2=s2(1,xc2);
    else
        st2=str2num(get(handles.edit_beg2,'String'));
    end
    
    % Create configuation file using output filename
    f=get(handles.edit_out,'String');
    fc=[f,'.conf'];
    fid=fopen(fc,'w');
    
    % Write signal filename(s) and related parameters
    if(get(handles.multi_sig,'Value')==1)
        sigs=handles.moresig;
        [m,n]=size(sigs);
        num=m/2+1;
        
        fprintf(fid,'%s',labels(12,:));
        for i=1:num
            if(i==1)
                sig=get(handles.edit_sig,'String');
            else
                sig=deblank(sigs(2*(i-1)-1,:));
            end
            fprintf(fid,' %s',sig);
        end
        fprintf(fid,'\n');
    else
        fprintf(fid,'%s %s\n',labels(12,:),f1); 
    end
    fprintf(fid,'%s %g\n',labels(6,:),st1);
    fprintf(fid,'%s %g\n',labels(7,:),e1);
    fprintf(fid,'%s %d\n\n',labels(8,:),nint1);
    
    % Write target filename(s) and related parameters
    if(get(handles.multi_sig,'Value')==1)
        sigs=handles.moresig;
        [m,n]=size(sigs);
        num=m/2+1;
        
        fprintf(fid,'%s',labels(13,:));
        for i=1:num
            if(i==1)
                targ=get(handles.edit_targ,'String');
            else
                targ=deblank(sigs(2*(i-1),:));
            end
            fprintf(fid,' %s',targ);
        end
        fprintf(fid,'\n');
    else
        fprintf(fid,'%s %s\n',labels(13,:),f2); 
    end
    fprintf(fid,'%s %g\n',labels(9,:),st2);
    fprintf(fid,'%s %g\n',labels(10,:),e2);
    fprintf(fid,'%s %d\n\n',labels(11,:),nint2);
    
    % Write penalty values
    fprintf(fid,'%s %d\n',labels(1,:),nmpen);
    fprintf(fid,'%s %d\n',labels(3,:),sppen);
    fprintf(fid,'targetspeed   %s\n',handles.target_speed);
    fprintf(fid,'%s %d\n',labels(2,:),spchpen);
    fprintf(fid,'%s %d\n',labels(5,:),tiepen);
    fprintf(fid,'%s %d\n\n',labels(4,:),gappen);
    
    
    fprintf(fid,'speeds %s\n\n',speeds2str(handles));
    
    % Write optional input files (tie and gap) if needed
    if(~isempty(handles.ties))
        ft=[f,'.tie'];
        if(xc1==1)
            t=handles.ties(:,[3 5]);
            save(ft,'t','-ASCII');
        else
            t=handles.ties(:,2:5);
            save(ft,'t','-ASCII');
        end
        fprintf(fid,'%s %s\n',labels(14,:),ft);
    end
    
    if(~isempty(handles.gaps))
        gaps=handles.gaps;
        for ax=1:2
            ind=find(gaps(:,2)==ax);
            if(~isempty(ind))
                if(ax==1)
                    fg=[str1 '_' f '.gap'];
                    fprintf(fid,'%s %s\n','series1gaps  ',fg);
                else
                    fg=[str2 '_' f,'.gap'];
                    fprintf(fid,'%s %s\n','series2gaps  ',fg);
                    s1=s2;
                    xc1=xc2;
                    yc1=yc2;
                end

                fgid=fopen(fg,'w');

                for i=1:size(gaps(ind,:))
                    left=gaps(ind(i),4);
                    right=gaps(ind(i),5);
                    a=num2str(left);
                    b=num2str(right);
                    if(yc1==2)  
                        core=0;
                        lend=min(s1(:,xc1));
                        rend=max(s1(:,xc1));
                    else
                        core=gaps(ind(i),3);
                        lend=min(s1(find(s1(:,1)==core),xc1));
                        rend=max(s1(find(s1(:,1)==core),xc1));
                    end
                    
                    if(left==lend)
                        a='NaN';
                    end
                    if(right==rend)
                        b='NaN';
                    end
                    fprintf(fgid,'%s %s %s\n',num2str(core), a, b);
                end
                fclose(fgid);
            end   
        end
    end
                
    % Write output filenames
    fm=[f,'.match'];
    fl=[f,'.log'];
    fprintf(fid,'%s %s\n',labels(15,:),fm);
    fprintf(fid,'%s %s\n',labels(16,:),fl);
    
    fclose(fid);
end

% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.
% --------------------------------------------------------------------
% Quit program (graphical user interface, match_gui)
ans=confirm(1);
if(ans==1)
    delete(handles.figure1);
end

% --------------------------------------------------------------------
function varargout = edit_beg1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_beg1.
% --------------------------------------------------------------------
% Set match begin point for signal

% Read in signal
f1=get(handles.edit_sig,'String');
n=0;
try
    s1=load(f1);
    [m n]=size(s1);
    xc=n-1;
    yc=n;
end

% If entered value not valid, set to beginning of signal
if(n==2|n==3)
    e=str2num(get(handles.edit_beg1,'String'));
    if(e<s1(1,xc) | e>s1(end,xc) | e>str2num(get(handles.edit_end1,'String')))
        e=min(s1(:,xc));
        set(handles.edit_beg1,'String',num2str(e));
        guidata(gcbo,handles);
    end
    axes(handles.axes1)
    hold on;
    delete(findobj('Marker','>','Parent',handles.axes1));
    x=find(s1(:,xc)<=e);
    plot(e,s1(x(end),yc),'g>','MarkerSize',7,'MarkerFaceColor','g');
    if(get(handles.pushbutton5,'Value')==1)
        set_axes(h,eventdata,handles,1);
    end
else
    disp('Warning: Invalid signal')
end

% --------------------------------------------------------------------
function varargout = edit_end1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_end1.
% --------------------------------------------------------------------
% Set match end point for signal

% Read in signal
f1=get(handles.edit_sig,'String');
n=0;
try
    s1=load(f1);
    [m n]=size(s1);
    xc=n-1;
    yc=n;
end

if(n==2|n==3)
    % If entered value not valid, set to end of signal
    e=str2num(get(handles.edit_end1,'String'));
    if(e<s1(1,xc) | e>s1(end,xc)| e<str2num(get(handles.edit_beg1,'String')))
        e=max(s1(:,xc));
        set(handles.edit_end1,'String',num2str(e));
        guidata(gcbo,handles);
    end
    axes(handles.axes1)
    hold on;
    delete(findobj('Marker','<','Parent',handles.axes1));
    x=find(s1(:,xc)<=e);
    plot(e,s1(x(end),yc),'g<','MarkerSize',7,'MarkerFaceColor','g');
    if(get(handles.pushbutton5,'Value')==1)
        set_axes(h,eventdata,handles,1);
    end
else
    disp('Warning: Invalid signal')
end

% --------------------------------------------------------------------
function varargout = edit_end2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_end2.
% --------------------------------------------------------------------
% Set match end point of target

% Read in target
f1=get(handles.edit_targ,'String');
n=0;
try
    s1=load(f1);
    [m n]=size(s1);
    xc=n-1;
    yc=n;
end

if(n==2|n==3)
    % If entered value not valid, set to end of target
    e=str2num(get(handles.edit_end2,'String'));
    if(e<s1(1,xc) | e>s1(end,xc)| e<str2num(get(handles.edit_beg2,'String')))
        e=max(s1(:,xc));
        set(handles.edit_end2,'String',num2str(e));
        guidata(gcbo,handles);
    end
    axes(handles.axes2)
    hold on;
    delete(findobj('Marker','<','Parent',handles.axes2));
    x=find(s1(:,xc)<=e);
    plot(e,s1(x(end),yc),'g<','MarkerSize',7,'MarkerFaceColor','g');
    if(get(handles.pushbutton5,'Value')==1)
        set_axes(h,eventdata,handles,1);
    end
else
    disp('Warning: Invalid target')
end

% --------------------------------------------------------------------
function varargout = edit_beg2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_beg2.
% --------------------------------------------------------------------
% Set match begin point of target

% Read in target
f1=get(handles.edit_targ,'String');
try
    s1=load(f1);
    [m n]=size(s1);
    xc=n-1;
    yc=n;
end

if(n==2|n==3)
    % If entered value not valid, set to beginning of target
    e=str2num(get(handles.edit_beg2,'String'));
    if(e<s1(1,xc) | e>s1(end,xc)| e>str2num(get(handles.edit_end2,'String')))
        e=min(s1(:,xc));
        set(handles.edit_beg2,'String',num2str(e));
        guidata(gcbo,handles);
    end
    axes(handles.axes2)
    hold on;
    delete(findobj('Marker','>','Parent',handles.axes2));
    x=find(s1(:,xc)<=e);
    plot(e,s1(x(end),yc),'g>','MarkerSize',7,'MarkerFaceColor','g');
    if(get(handles.pushbutton5,'Value')==1)
        set_axes(h,eventdata,handles,1);
    end
else
    disp('Warning: Invalid target')
end

% --------------------------------------------------------------------
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton5.
% --------------------------------------------------------------------
% Button to set axes of signal and target to the match begin and end points
set_axes(h,eventdata,handles,get(handles.pushbutton5,'Value'));

% --------------------------------------------------------------------
function set_axes(h, eventdata, handles, value)
% --------------------------------------------------------------------
% Set axes of signal and target to the match begin and end points

handles=guidata(handles.figure1);
% Read in signal and target
f1=get(handles.edit_sig,'String');
try
    s1=load(f1);
    [m n]=size(s1);
    xc1=n-1;
end

f1=get(handles.edit_targ,'String');
try
    s2=load(f1);
    [m n]=size(s2);
    xc2=n-1;
end

if(value==1)
    % If match begin and end values already entered, set axes
    axes(handles.axes1)
    st1=str2num(get(handles.edit_beg1,'String'));
    e1=str2num(get(handles.edit_end1,'String'));
    sub1=s1(find(s1(:,xc1)>=st1 & s1(:,xc1)<=e1),:);
    xlim([st1 e1])
    ylim(handles.axes1,'auto');
    axes(handles.axes2)
    st2=str2num(get(handles.edit_beg2,'String'));
    e2=str2num(get(handles.edit_end2,'String'));
    sub2=s2(find(s2(:,xc2)>=st2 & s2(:,xc2)<=e2),:);
    xlim([st2 e2])
    ylim(handles.axes1,'auto');
    
else
    % If match begin and end values don't exist, set to extrema of 
    % signal and target and set axes
    axes(handles.axes1)
    f1=get(handles.edit_sig,'String');
    s1=load(f1);
    st1=min(s1(:,xc1));
    e1=max(s1(:,xc1));
    xlim([st1 e1])
    ylim(handles.axes1,'auto');
    axes(handles.axes2)
    f1=get(handles.edit_targ,'String');
    s2=load(f1);
    st2=min(s2(:,xc2));
    e2=max(s2(:,xc2));
    xlim([st2 e2])
    ylim(handles.axes2,'auto');
end

% set slider bars values
set(handles.slider1,'min',min(s1(:,xc1)),'max',max(s1(:,xc1)),'Value',st1,'sliderstep',[0 0]);
set(handles.slider2,'min',min(s2(:,xc2)),'max',max(s2(:,xc2)),'Value',st2,'sliderstep',[0 0]);
% put axes limits in text boxes
axes(handles.axes1)
b=axis;
c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
set(handles.edit_axeslim1,'String',c);
axes(handles.axes2)
b=axis;
c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
set(handles.edit_axeslim2,'String',c);

% --------------------------------------------------------------------
function varargout = push_zoom_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_zoom.
% --------------------------------------------------------------------
% Reduce x-axis limits of signal to zoom in
% The left-hand value is constant, the total axis size is reduced by 50%

x=xlim(handles.axes1);
xlim(handles.axes1, [x(1) x(1)+(x(2)-x(1))/2]);
ylim(handles.axes1,'auto');
set(handles.slider1,'Value',x(1),'SliderStep',[.03 .06]);
axes(handles.axes1)
b=axis;
c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
set(handles.edit_axeslim1,'String',c);

% If synchronize button is depressed, set target (and results) axes to be
% the same as the signal's
if(get(handles.transfer_axeslim,'Value')==1)
    xlim(handles.axes2, [b(1) b(2)]);
    ylim(handles.axes2, [b(3) b(4)]);
    if(b(1)>get(handles.slider2,'Max'))
        set(handles.slider2,'Value',get(handles.slider2,'Max'));
    elseif(b(1)<get(handles.slider2,'Min'))    
        set(handles.slider2,'Value',get(handles.slider2,'Min'));
    else
        set(handles.slider2,'Value',b(1),'SliderStep',[.03 .06]);
    end
    set(handles.edit_axeslim2,'String',c);
    if(handles.match_fig~=0)
        h=handles.results_axes;
        for i=1:max(size(h))
            xlim(handles.results_axes(i), [b(1) b(2)]);
        end
        if(get(handles.multi_sig,'Value')==0)
            ylim(handles.results_axes(1), [b(3) b(4)]);
        end
        set(handles.results_lim,'String',c);
    end
end


% --------------------------------------------------------------------
function varargout = push_zoom2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_zoom.
% --------------------------------------------------------------------
% Reduce x-axis limits of target to zoom in
% The left-hand value is constant, the total axis size is reduced by 50%

x=xlim(handles.axes2);
xlim(handles.axes2, [x(1) x(1)+(x(2)-x(1))/2]);
ylim(handles.axes2,'auto');
set(handles.slider2,'Value',x(1),'SliderStep',[.03 .06]);
axes(handles.axes2)
b=axis;
c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
set(handles.edit_axeslim2,'String',c);

% If synchronize button is depressed, set signal (and results) axes to be
% the same as the target's
if(get(handles.transfer_axeslim,'Value')==1)
    xlim(handles.axes1, [b(1) b(2)]);
    ylim(handles.axes1, [b(3) b(4)]);
    if(b(1)>get(handles.slider1,'Max'))
        set(handles.slider1,'Value',get(handles.slider1,'Max'));
    elseif(b(1)<get(handles.slider1,'Min'))    
        set(handles.slider1,'Value',get(handles.slider1,'Min'));
    else
        set(handles.slider1,'Value',b(1),'SliderStep',[.03 .06]);
    end
    set(handles.edit_axeslim1,'String',c);
    if(handles.match_fig~=0)
        h=handles.results_axes;
        for i=1:max(size(h))
            xlim(handles.results_axes(i), [b(1) b(2)]);
        end
        if(get(handles.multi_sig,'Value')==0)
            ylim(handles.results_axes(1), [b(3) b(4)]);
        end
        set(handles.results_lim,'String',c);
    end
end


% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
% --------------------------------------------------------------------
% Slider setting is used to adjust the axes of signal plot.

x=xlim(handles.axes1);
f1=get(handles.edit_sig,'String');
s1=load(f1);
x1max=max(s1(:,end-1));
range=x1max-min(s1(:,end-1));
if(x(2)-x(1) > range)
    range=x(2)-x(1);
end
set(handles.slider1,'Min',min(s1(:,end-1)),'Max',x1max,'sliderstep',[.03 .06]);

%if(x(1) < x1max)
    v=get(handles.slider1,'Value')
    xlim(handles.axes1, [v v+x(2)-x(1)]);
    ylim(handles.axes1,'auto');
    axes(handles.axes1)
    b=axis;
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim1,'String',c);

% If synchronize button is depressed, set target (and results) axes to be
% the same as the signal's
    if(get(handles.transfer_axeslim,'Value')==1)
        xlim(handles.axes2, [b(1) b(2)]);
        ylim(handles.axes2, [b(3) b(4)]);
        if(b(1)>get(handles.slider2,'Max'))
            set(handles.slider2,'Value',get(handles.slider2,'Max'));
        elseif(b(1)<get(handles.slider2,'Min'))    
            set(handles.slider2,'Value',get(handles.slider2,'Min'));
        else
            set(handles.slider2,'Value',b(1),'SliderStep',[.03 .06]);
        end
        set(handles.edit_axeslim2,'String',c);
        if(handles.match_fig~=0)
            h=handles.results_axes;
            for i=1:max(size(h))
                xlim(handles.results_axes(i), [b(1) b(2)]);
            end
            ylim(handles.results_axes(1), [b(3) b(4)]);
            set(handles.results_lim,'String',c);
        end
    end
    %end

% --------------------------------------------------------------------
function varargout = slider2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider2.
% --------------------------------------------------------------------
% Slider setting is used to adjust the axes of target plot.

x=xlim(handles.axes2);
f1=get(handles.edit_targ,'String');
s1=load(f1);
x2max=max(s1(:,end-1));
range=x2max-min(s1(:,end-1));
if(x(2)-x(1) > range)
    range=x(2)-x(1);
end
set(handles.slider2,'Min',min(s1(:,end-1)),'Max',x2max,'sliderstep',[.03 .06]);

%if(x(1) < x2max)
    v=get(handles.slider2,'Value');
    xlim(handles.axes2, [v v+x(2)-x(1)]);
    ylim(handles.axes2,'auto');
    axes(handles.axes2)
    b=axis;
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim2,'String',c);
    
% If synchronize button is depressed, set signal (and results) axes to be
% the same as the target's
    if(get(handles.transfer_axeslim,'Value')==1)
        xlim(handles.axes1, [b(1) b(2)]);
        ylim(handles.axes1, [b(3) b(4)]);
        if(b(1)>get(handles.slider1,'Max'))
            set(handles.slider1,'Value',get(handles.slider1,'Max'));
        elseif(b(1)<get(handles.slider1,'Min'))    
            set(handles.slider1,'Value',get(handles.slider1,'Min'));
        else
            set(handles.slider1,'Value',b(1),'SliderStep',[.03 .06]);
        end
        set(handles.edit_axeslim1,'String',c);
        if(handles.match_fig~=0)
            h=handles.results_axes;
            for i=1:max(size(h))
                xlim(handles.results_axes(i), [b(1) b(2)]);
            end
            ylim(handles.results_axes(1), [b(3) b(4)]);
            set(handles.results_lim,'String',c);
        end
    end
    %end
 
% --------------------------------------------------------------------
function varargout = radio_tie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_tie.
% --------------------------------------------------------------------
% If this radio button is selected, buttons above plots control tie points

if(get(handles.radio_tie,'Value')==1)
    set(handles.radio_gap,'Value',0);
else
    set(handles.radio_gap,'Value',1);
end

% --------------------------------------------------------------------
function varargout = radio_gap_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_gap.
% --------------------------------------------------------------------
% If this radio button is selected, buttons above plots control inserted gaps

if(get(handles.radio_gap,'Value')==1)
    set(handles.radio_tie,'Value',0);
else
    set(handles.radio_tie,'Value',1);
end

%--------------------------------------------------------------------
function varargout = toggle_tie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.toggle_tie.
%--------------------------------------------------------------------
% This button prompts the user to add new tie points or gaps 
% until the 'Enter' key is pressed

% Read in signal and target
handles=guidata(handles.figure1);
f1=get(handles.edit_sig,'String');
f2=get(handles.edit_targ,'String');
try
    s1=load(f1);
    [m n]=size(s1);
    xc1=n-1;
    yc1=n;
end
try
    s2=load(f2);
    [m n]=size(s2);
    xc2=n-1;
    yc2=n;
end

% Test for valid number of columns in signal and target files
if((xc1==1 | xc1==2) & (xc2==1 | xc2==2))
    
    if(get(handles.toggle_tie,'Value')==1)
        core1=0;
        core2=0;
        set(handles.toggle_tie,'String','Press Enter to Stop');    
        
        % Add tie point
        if(get(handles.radio_tie,'Value')==1 & get(handles.toggle_tie,'Value')==1)
            % Get x-value (depth/time) of tie point in signal
            [x1,y1]=ginput(1);
            % Continue only if x-value is valid
            if(~isempty(x1) & x1>=min(s1(:,xc1)) & x1<=max(s1(:,xc1))& get(gcf,'CurrentAxes')==handles.axes1)
                axes(handles.axes1)
                val=find(s1(:,xc1)<=x1);
                % If multi-core signal, test if point falls on overlapping cores
                if(yc1>2)
                    a=size(val);
                    b=size([1:val(end)]');
                    c=s1(val(end),1);
                    if(a(1)~=b(1))
                        c1=min(s1(find(s1(:,1)==c),xc1));
                        c2=max(s1(find(s1(:,1)==c),xc1));
                        ind=find(s1(val,1)~=c);
                        d=s1(val(ind(end)),1);
                        d1=min(s1(find(s1(:,1)==d),xc1));
                        d2=max(s1(find(s1(:,1)==d),xc1));
                        v=[c c1 c2 d d1 d2];
                        % ask user in which core to place tie point
                        ans=pick_core(v);
                        if(ans==1)
                            core1=d;
                            val=val(ind);
                        else
                            core1=c;
                        end
                    else
                        core1=c;
                    end
                end
                
                % Plot and label tie point in signal
                plot(x1,s1(val(end),yc1),'k+','MarkerSize',5);
                text(x1,s1(val(end),yc1),num2str(handles.nties),'VerticalAlignment','bottom',...
                    'HorizontalAlignment','center','Color','k','FontSize',12);
                
                % Prompt for corresponding point in target
                [x2,y2]=ginput(1);
                
                % Continue if point is in the target, otherwise delete marker in signal
                if(~isempty(x2) & x2>=min(s2(:,xc2)) & x2<=max(s2(:,xc2)) & ...
                        get(gcf,'CurrentAxes')==handles.axes2)
                    axes(handles.axes2)
                    val=find(s2(:,xc2)<=x2);
                    % If multi-core signal, test if point falls on overlapping cores
                    if(yc2>2)
                        a=size(val);
                        b=size([1:val(end)]');
                        c=s2(val(end),1);
                        if(a(1)~=b(1))
                            c1=min(s2(find(s2(:,1)==c),xc2));
                            c2=max(s2(find(s2(:,1)==c),xc2));
                            ind=find(s2(val,1)~=c);
                            d=s2(val(ind(end)),1);
                            d1=min(s2(find(s2(:,1)==d),xc2));
                            d2=max(s2(find(s2(:,1)==d),xc2));
                            v=[c c1 c2 d d1 d2];
                            ans=pick_core(v);
                            if(ans==1)
                                core2=d;
                                val=val(ind);
                            else
                                core2=c;
                            end
                        else
                            core2=c;
                        end 
                    end
                    
                    % Plot tie point in target
                    if(get(handles.radio_tie,'Value')==1)
                        plot(x2,s2(val(end),yc2),'k+','MarkerSize',5);
                        text(x2,s2(val(end),yc2),num2str(handles.nties),'VerticalAlignment','bottom',...
                            'HorizontalAlignment','center','Color','k','FontSize',12);
                        % Save tie point
                        handles.ties(end+1,:)=[handles.nties core1 x1 core2 x2];
                        handles.nties=handles.nties+1;    
                    end            
                else
                    % Deleting tie point because invalid
                    set(handles.toggle_tie,'Value',0);
                    set(handles.toggle_tie,'String','Add');    
                    hobj=findobj('String',num2str(handles.nties),'Parent',handles.axes1,'Color','k');
                    x=get(hobj,'Position');
                    delete(hobj);
                    delete(findobj('XData',x(1),'Parent',handles.axes1));
                end;
            else
                % If enter key pressed instead of clicking on graph, stop adding gaps  
                set(handles.toggle_tie,'Value',0);
                set(handles.toggle_tie,'String','Add');    
            end            
            guidata(gcbo,handles);
            
        % If gap radio button selected, add a gap to signal or target 
        elseif(get(handles.toggle_tie,'Value')==1 & get(handles.radio_tie,'Value')==0)
            guidata(gcbo,handles);
            % Prompt for gap position
            [x1,y1]=ginput(1);
            % Get gap description and plot gap
            gap=plot_gaps(h, eventdata, handles, -1, x1);
            if(gap(1)~=0)
                % If valid, save gap
                handles.gaps(end+1,:)=[handles.ngaps gap];
                handles.ngaps=handles.ngaps+1;
                guidata(gcbo,handles);
            else
                x1=[];
            end
        end
    end
    
    % If enter key pressed instead of clicking on graph, stop adding gaps  
    if(isempty(x1))
        set(handles.toggle_tie,'Value',0);  
        set(handles.toggle_tie,'String','Add');        
    end
    
    % Prompt for a new tie point/gap
    if(get(handles.toggle_tie,'Value')==1)
        toggle_tie_Callback(h, eventdata, handles, varargin);
    end
else
    % Change button label back to 'Add'
    set(handles.toggle_tie,'String','Add');    
end

% --------------------------------------------------------------------
function varargout = load_tie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.load_tie.
% --------------------------------------------------------------------
% Load tie points/gaps from a file (deleting any previously existing tie points/gaps)

if(get(handles.radio_tie,'Value')==1)
 % Load tie points
    tiefile=getconf(2);
    if(~strcmp(tiefile,'cancel'))
        t=load(tiefile);
        [m n]=size(t);
        ind=[1:m]';
        % Save tie points in memory
        if(n==2)
            handles.ties=[ind 0*ind t(:,1) 0*ind t(:,2)];
        else
            handles.ties=[ind t];
        end
        handles.nties=m+1;
        guidata(gcbo,handles);
        plot_ties(h, eventdata, handles,varargin);
    end
else   % Load gaps
    ax=confirm(4)+1;     % Gaps for signal or target? 1=signal 2=target
    gapfile=getconf(3);  % Get filename
    g=[];
    g=load(gapfile);
    if(~strcmp(gapfile,'cancel'))
        if(ax==1)
            f1=get(handles.edit_sig,'String');
            delete(findobj('Parent',handles.axes1,'Color','r'));
            handles.gaps(find(handles.gaps(:,2)==1),:)=[];
        else
            f1=get(handles.edit_targ,'String');
            delete(findobj('Parent',handles.axes2,'Color','r'));
            handles.gaps(find(handles.gaps(:,2)==2),:)=[];
        end
        if(~isempty(g))
            try
                s1=load(f1);
                [m n1]=size(s1);
            end
            [m n]=size(g);
            for i=1:m    % Loop over all gaps
                x=g(i,2:3);   
                % Get start and end points of gaps
                if(n1==2)   % For series without cores
                    if(isnan(x(1)))
                        x(1)=min(s1(:,1));
                    end
                    if(isnan(x(2)))
                        x(2)=max(s1(:,1));
                    end
                else        % For series with cores
                    if(isnan(x(1)))
                        x(1)=min(s1(find(s1(:,1)==g(i,1)),n1-1));
                    end
                    if(isnan(x(2)))
                        x(2)=max(s1(find(s1(:,1)==g(i,1)),n1-1));
                    end
                end
                % Save gaps in memory
                handles.gaps(end+1,:)=[handles.ngaps ax g(i,1) x(1) x(2)];
                handles.ngaps=handles.ngaps+1;
            end
        end
        guidata(gcbo,handles);
        plot_ties(h, eventdata, handles, 1);   % Plot gaps
    end
end    

% --------------------------------------------------------------------
function answer=plot_gaps(h, eventdata, handles, num, pos)
% --------------------------------------------------------------------
% Add or modify (and plot) inserted gap
% integer num is the gap label if pre-existing gap or -1 if the gap is new
% pos tells the type of gap/ how it should be moved
% 1=core top/start, 2=core bottom/end, 3=core middle (7 or 8), 4=single point
% 5=cancel,         6=entire core,     7=move gap top,         8=move gap bottom
% answer is an array identiying the series#, core#, start, and end of gap
% Gaps are plotted as red lines bounded by red squares labeled with red numbers

cancel=0;

% Load signal or target
if(get(gcf,'CurrentAxes')==handles.axes1)
    f1=get(handles.edit_sig,'String');
    try
        s1=load(f1);
        [m n]=size(s1);
        xc1=n-1;
        yc1=n;
        ax=1;
    end
else
    f2=get(handles.edit_targ,'String');
    try
        s1=load(f2);
        [m n]=size(s1);
        xc1=n-1;
        yc1=n;
        ax=2;
    end
end

% Get gap position from argument list (adding gap) or user (moving gap)
if(num==-1)
    x1=pos;
    g=[];
else
    [x1,y1]=ginput(1);
    g=find(handles.gaps(:,1)==num);
end

% If gap is valid ...
if(~isempty(x1)& x1>=min(s1(:,xc1)) & x1<=max(s1(:,xc1)))
    val=find(s1(:,xc1)<=x1);
    if(yc1==2)
        core1=0;
    else
        % Test for overlapping cores
        a=size(val);
        b=size([1:val(end)]');
        c=s1(val(end),1);
        if(a(1)~=b(1))
            c1=min(s1(find(s1(:,1)==c),xc1));
            c2=max(s1(find(s1(:,1)==c),xc1));
            ind=find(s1(val,1)~=c);
            d=s1(val(ind(end)),1);
            d1=min(s1(find(s1(:,1)==d),xc1));
            d2=max(s1(find(s1(:,1)==d),xc1));
            v=[c c1 c2 d d1 d2];
            ans=pick_core(v);
            if(ans==1)
                core1=d;
                val=val(ind);
            else
                core1=c;
            end
        else
            core1=c;
        end
    end
    
    if(isempty(g))
        % For new gaps, ask for type of gap and get new label number
        ans=gap_size;
        label=num2str(handles.ngaps);
    else
        % For moving gaps, change appropriate property of gap
        if(pos==7)                %moving front end of gap
            x2=handles.gaps(g,5);
        elseif (pos==8)           %moving back end of gap
            x2=handles.gaps(g,4);
            pos=7;
        else
            % current gap position
            left=handles.gaps(g,4);
            right=handles.gaps(g,5);
            % find valid gap limits (i.e. core/series ends
            if(yc1==2)
                lend=min(s1(:,xc1));
                rend=max(s1(:,xc1));
            else
                core=handles.gaps(g,3);
                lend=min(s1(find(s1(:,1)==core),xc1));
                rend=max(s1(find(s1(:,1)==core),xc1));
            end
            % identify type of gap and set pos
            if(left==right)
                pos=4
            elseif(left==lend)
                if(right==rend)
                    pos=6
                else
                    pos=1
                end
            elseif(right==rend)
                pos=2
            end
        end
        ans=pos;
        label=num2str(num);
    end
    
    % Set appropriate start and end values for gap (plot end if necessary)
    % x1=clicked position/gap start, x2=gap end
    % val=nearest preceding data point (for plotting purposes)
    switch ans
    case 1                      % gap=[core top, x1]
        x2=x1;
        if(yc1==2)
            x1=min(s1(:,xc1));
        else
            x1=min(s1(find(s1(:,1)==core1),xc1));
        end
        val=find(s1(:,xc1)<=x1);
    case 2                      % gap=[x1, core end]
        if(yc1==2)
            x2=max(s1(:,xc1));
        else
            x2=max(s1(find(s1(:,1)==core1),xc1));
        end
    case 3                      % gap=[x1,x2] or [x2,x1]
        plot(x1,s1(val(end),yc1),'rs','MarkerSize',4);
        [x2,y2]=ginput(1);
        if(~isempty(x2))
            if(yc1==2)
                if(x2<min(s1(:,xc1)))
                    x2=min(s1(:,xc1));
                elseif (x2>max(s1(:,xc1)))
                    x2=max(s1(:,xc1));                                  
                else
                    x2=s1(max(find(s1(:,xc1)<=x2)),xc1);
                end
            else
                ind=find(s1(:,1)==core1);
                if(x2<min(s1(ind,xc1)))
                    x2=min(s1(ind,xc1));
                elseif (x2>max(s1(ind,xc1)))
                    x2=max(s1(ind,xc1));                                 
                else
                    i=find(s1(:,xc1)<=x2 & s1(:,1)==core1);
                    x2=s1(i(end),xc1);
                end
            end
        else
            cancel=1;
        end
    case 4                       % point gap x1=x2
        x2=x1;
    case 6                       % gap=[core top,core end]
        if(yc1==2)
            x1=min(s1(:,xc1));
            x2=max(s1(:,xc1));
        else
            x1=min(s1(find(s1(:,1)==core1),xc1));
            x2=max(s1(find(s1(:,1)==core1),xc1));
        end
        val=find(s1(:,xc1)<=x1);
    case 7                       % gap=[x1,x2=old gap end]
        % Delete old markers of start and end
        oldx=get(findobj('String',num2str(num),'Parent',handles.axes1,'Color','r'),'Position');
        delete(findobj('XData',oldx(1),'Parent',handles.axes1,'Color','r'));
        delete(findobj('String',num2str(num),'Parent',handles.axes1,'Color','r'));
        
        plot(x1,s1(val(end),yc1),'rs','MarkerSize',4);
        if(yc1==2)
            if(x2<min(s1(:,xc1)))
                x2=min(s1(:,xc1));
            elseif (x2>max(s1(:,xc1)))
                x2=max(s1(:,xc1));                                  
            end
        else
            ind=find(s1(:,1)==core1);
            if(x2<min(s1(ind,xc1)))
                x2=min(s1(ind,xc1));
            elseif (x2>max(s1(ind,xc1)))
                x2=max(s1(ind,xc1));                                 
            end
        end
        
    otherwise
        cancel=1;
    end
    
    % Plot gap and return gap info (answer)
    if(~cancel)
        if(x2<x1)       %Flip x1 and x2 if x1>x2
            a=x1;
            x1=x2;
            x2=a;
            val=find(s1(:,xc1)<=x1);
        end
        if(ans~=1)
            plot(x1,s1(val(end),yc1),'rs','MarkerSize',4);
            text(x1,s1(val(end),yc1),label,'VerticalAlignment','bottom',...
                'HorizontalAlignment','center','Color','r','FontSize',12);
        end
        if(yc1==2)
            val2=find(s1(:,xc1)<=x2);
        else
            val2=find(s1(:,1)==core1 & s1(:,xc1)<=x2);
        end
        if(ans==1 | ans==3 | ans==7)
            plot(x2,s1(val2(end),yc1),'rs','MarkerSize',4);
            text(x2,s1(val2(end),yc1),label,'VerticalAlignment','bottom',...
                'HorizontalAlignment','center','Color','r','FontSize',12);
        end
        if(ans~=4)
            plot(s1(val(end):val2(end),xc1),s1(val(end):val2(end),yc1),'r');
        end
        answer=[ax core1 x1 x2];
    else
        answer=[0 0 0 0];
    end
else
        answer=[0 0 0 0];
end


% --------------------------------------------------------------------
function varargout = del_tie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.del_tie.
% --------------------------------------------------------------------
% Remove a tie point or gap from graphs and memory

% Get number of tie point or gap from text box
nt=get(handles.del_tie,'String');
if(~isempty(nt))
    if(get(handles.radio_tie,'Value')==1)
        % Find tie point labeled by nt
        if(isempty(handles.ties))
            disp(['Tie Point ' nt ' does not exist'])
            set(handles.del_tie,'String','');
            return
        end
        ind=find(handles.ties(:,1)==str2num(nt));
        if(~isempty(ind))
            % Delete tie point nt from memory
            handles.ties(ind,:)=[];
            guidata(gcbo,handles);
            % Remove tie point nt from graph
            hobj1=findobj('String',nt,'Parent',handles.axes1,'Color','k');
            hobj2=findobj('String',nt,'Parent',handles.axes2,'Color','k');
            x=get(hobj1,'Position');
            delete(hobj1);
            delete(findobj('XData',x(1),'Parent',handles.axes1));
            y=get(hobj2,'Position');
            delete(findobj('XData',y(1),'Parent',handles.axes2));
            delete(hobj2);
        else
            disp(['Tie point ' nt ' does not exist'])
            set(handles.del_tie,'String','');
        end
    else
        % Find gap labeled by nt
        if(isempty(handles.gaps))
            disp(['Gap ' nt ' does not exist'])
            set(handles.del_tie,'String','');
            return
        end
        ind=find(handles.gaps(:,1)==str2num(nt));
        if(isempty(ind))
            disp(['Gap ' nt ' does not exist'])
            set(handles.del_tie,'String','');
        else
            %Get series/axes (signal or target), start, and end of gap
            ax=handles.gaps(ind,2);
            x1=handles.gaps(ind,4);
            x2=handles.gaps(ind,5);
            % Delete gap from memory
            handles.gaps(ind,:)=[];
            guidata(gcbo,handles);
            % Erase gap from plot
            if(ax==1)               % gap in signal
                f1=get(handles.edit_sig,'String');
                try
                    s1=load(f1);
                    [m n]=size(s1);
                    xc1=n-1;
                    yc1=n;
                end
                axes(handles.axes1)
                val=find(s1(:,xc1)<=x1);
                val2=find(s1(:,xc1)<=x2);
                hobj1=findobj('String',nt,'Parent',handles.axes1,'Color','r');
                delete(hobj1);
                delete(findobj('XData',x1,'Parent',handles.axes1));
                delete(findobj('XData',x2,'Parent',handles.axes1));
                plot(s1(val(end):val2(end),xc1),s1(val(end):val2(end),yc1),'m');
                
            else                        % gap in target
                f1=get(handles.edit_targ,'String');
                try
                    s1=load(f1);
                    [m n]=size(s1);
                    xc1=n-1;
                    yc1=n;
                end
                axes(handles.axes2)
                val=find(s1(:,xc1)<=x1);
                val2=find(s1(:,xc1)<=x2);
                hobj2=findobj('String',nt,'Parent',handles.axes2,'Color','r');
                delete(findobj('XData',x1,'Parent',handles.axes2));
                delete(findobj('XData',x2,'Parent',handles.axes2));
                delete(hobj2);
                plot(s1(val(end):val2(end),xc1),s1(val(end):val2(end),yc1),'Color',[.35 .35 1]);
            end
        end
    end
    set(handles.del_tie,'String','');
end

% --------------------------------------------------------------------
function varargout = button_move_tie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.button_move_tie.
% --------------------------------------------------------------------
% When button is depressed, clicking on tie points and gaps will 
% allow the user to move them

if(get(handles.button_move_tie,'Value')==0)
    set(handles.button_move_tie,'String','Move');   
else
    set(handles.button_move_tie,'String','Press to lock');   
end

% --------------------------------------------------------------------
function varargout = text1_ButtondownFcn(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Detects clicks on tie point/gap labels in signal and prompts user to move them

f1=get(handles.edit_sig,'String');
try
    s1=load(f1);
    [m n]=size(s1);
    xc1=n-1;
    yc1=n;
end
num=get(gcbo,'String');   % num = tie point (or gap) label number
pos=get(gcbo,'Position'); % pos = position of tie point (gap)/label

% Only moves point if Move button is depressed
if(get(handles.button_move_tie,'Value')==1)
    % If tie points are selected and label was a tie point label (i.e. black)
    if(get(handles.radio_tie,'Value')==1 & get(gcbo,'Color')==[0 0 0])
        [x1,y1]=ginput(1);         %Get new position
        if(~isempty(x1) & x1>=min(s1(:,xc1)) & x1<=max(s1(:,xc1)))
            axes(handles.axes1)
            % Remove old tie point from graph
            ind=find(handles.ties(:,1)==str2num(num));
            delete(findobj('String',num,'Parent',handles.axes1,'Color','k'));
            delete(findobj('XData',pos(1),'Parent',handles.axes1));
            val=find(s1(:,xc1)<=x1);
            core1=handles.ties(ind,2);
            % check for overlapping cores at new tie point position 
            if(yc1>2)
                a=size(val);
                b=size([1:val(end)]');
                c=s1(val(end),1);
                if(a(1)~=b(1))
                    c1=min(s1(find(s1(:,1)==c),xc1));
                    c2=max(s1(find(s1(:,1)==c),xc1));
                    ind1=find(s1(val,1)~=c);
                    d=s1(val(ind1(end)),1);
                    d1=min(s1(find(s1(:,1)==d),xc1));
                    d2=max(s1(find(s1(:,1)==d),xc1));
                    v=[c c1 c2 d d1 d2];
                    ans=pick_core(v);
                    if(ans==1)
                        core1=d;
                        val=val(ind1);
                    else
                        core1=c;
                    end
                else
                    core1=c;
                end
            end
            % Store and plot new tie point
            handles.ties(ind,2)=core1;
            handles.ties(ind,3)=s1(val(end),xc1);
            plot(x1,s1(val(end),yc1),'k+','MarkerSize',5);
            text(x1,s1(val(end),yc1),num,'VerticalAlignment','bottom',...
                'HorizontalAlignment','center','Color','k','FontSize',12);
        end
        
    % If gap is selected and label is a gap label (red)    
    elseif(get(handles.radio_gap,'Value')==1 & get(gcbo,'Color')==[1 0 0])
        n=str2num(num);                   % gap label number
        ind=find(handles.gaps(:,1)==n);   % gap index
        hobj1=findobj('String',num,'Parent',handles.axes1,'Color','r');
        s=size(hobj1);     
        if(s(1)>1)          % if multiple labels (i.e. gap in core middle)
            % remove gap from graph
            x1=get(hobj1(1),'Position');
            x2=get(hobj1(2),'Position');
            delete(findobj('String',num,'Parent',handles.axes1,'Color','r','Position',pos));
            delete(findobj('XData',pos(1),'Parent',handles.axes1));
            val=find(s1(:,xc1)<=handles.gaps(ind,4));
            val2=find(s1(:,xc1)<=handles.gaps(ind,5));
            plot(s1(val(end):val2(end),xc1),s1(val(end):val2(end),yc1),'m');
            % get new gap position & plot it
            if(min([x1(1) x2(1)])==pos(1))
                % move gap start
                v=plot_gaps(h, eventdata, handles, n, 7);
            else
                % move gap end
                v=plot_gaps(h, eventdata, handles, n, 8);
            end
        else               % if gap has only one moveable end
            % remove gap from graph
            delete(findobj('String',num,'Parent',handles.axes1,'Color','r'));
            delete(findobj('XData',pos(1),'Parent',handles.axes1));
            val=find(s1(:,xc1)<=handles.gaps(ind,4));
            val2=find(s1(:,xc1)<=handles.gaps(ind,5));
            plot(s1(val(end):val2(end),xc1),s1(val(end):val2(end),yc1),'m');
            % get new gap position & plot it
            v=plot_gaps(h, eventdata, handles, n, 1);
        end
        %store modified gap in memory    
        handles.gaps(ind,:)=[n v];
    end
end
guidata(handles.figure1,handles);

% --------------------------------------------------------------------
function varargout = text2_ButtondownFcn(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Detects clicks on tie point/gap labels in target and prompts user to move them
% see comments in code for function text1_ButtondownFcn

f1=get(handles.edit_targ,'String');
try
    s1=load(f1);
    [m n]=size(s1);
    xc2=n-1;
    yc2=n;
end
num=get(gcbo,'String');
pos=get(gcbo,'Position');

if(get(handles.button_move_tie,'Value')==1)
    if(get(handles.radio_tie,'Value')==1 & get(gcbo,'Color')==[0 0 0])
        [x1,y1]=ginput(1);
        if(~isempty(x1) & x1>=min(s1(:,xc2)) & x1<=max(s1(:,xc2)))
            axes(handles.axes2)
            ind=find(handles.ties(:,1)==str2num(num));
            delete(findobj('String',num,'Parent',handles.axes2,'Color','k'));
            delete(findobj('XData',pos(1),'Parent',handles.axes2));
            val=find(s1(:,xc2)<=x1);
            core2=handles.ties(ind,4);
            if(yc2>2)
                a=size(val);
                b=size([1:val(end)]');
                c=s1(val(end),1);
                if(a(1)~=b(1))
                    c1=min(s1(find(s1(:,1)==c),xc2));
                    c2=max(s1(find(s1(:,1)==c),xc2));
                    ind1=find(s1(val,1)~=c);
                    d=s1(val(ind1(end)),1);
                    d1=min(s1(find(s1(:,1)==d),xc2));
                    d2=max(s1(find(s1(:,1)==d),xc2));
                    v=[c c1 c2 d d1 d2];
                    ans=pick_core(v);
                    if(ans==1)
                        core2=d;
                        val=val(ind1);
                    else
                        core2=c;
                    end
                else
                    core2=c;
                end
            end
            handles.ties(ind,4)=core2;
            handles.ties(ind,5)=s1(val(end),xc2);
            plot(x1,s1(val(end),yc2),'k+','MarkerSize',5);
            text(x1,s1(val(end),yc2),num,'VerticalAlignment','bottom',...            
                'HorizontalAlignment','center','Color','k','FontSize',12);
        end
        
    elseif(get(handles.radio_gap,'Value')==1 & get(gcbo,'Color')==[1 0 0])
        n=str2num(num);
        ind=find(handles.gaps(:,1)==n);
        hobj1=findobj('String',num,'Parent',handles.axes2,'Color','r');
        s=size(hobj1);
        if(s(1)>1)
            x1=get(hobj1(1),'Position');
            x2=get(hobj1(2),'Position');
            delete(findobj('String',num,'Parent',handles.axes2,'Color','r','Position',pos));
            delete(findobj('XData',pos(1),'Parent',handles.axes2));
            val=find(s1(:,xc2)<=handles.gaps(ind,4));
            val2=find(s1(:,xc2)<=handles.gaps(ind,5));
            plot(s1(val(end):val2(end),xc2),s1(val(end):val2(end),yc2),'Color',[.35 .35 1]);
            if(min([x1(1) x2(1)])==pos(1))
                v=plot_gaps(h, eventdata, handles, n, 7);
            else
                v=plot_gaps(h, eventdata, handles, n, 8);
            end
        else
            delete(findobj('String',num,'Parent',handles.axes2,'Color','r'));
            delete(findobj('XData',pos(1),'Parent',handles.axes2));
            val=find(s1(:,xc2)<=handles.gaps(ind,4));
            val2=find(s1(:,xc2)<=handles.gaps(ind,5));
            plot(s1(val(end):val2(end),xc2),s1(val(end):val2(end),yc2),'Color',[.35 .35 1]);
            v=plot_gaps(h, eventdata, handles, n, 1);
        end
            
        handles.gaps(ind,:)=[n v];
    end
end
guidata(handles.figure1,handles);

% --------------------------------------------------------------------
function varargout = push_refresh_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_refresh.
% --------------------------------------------------------------------
% Renumbers all tie points or gaps to be numbered sequentially
% Also displays a warning if tie points are crossed 

% If tie points selected
if(get(handles.radio_tie,'Value')==1)
    [m,n]=size(handles.ties);
    handles.nties=m+1;
    % Delete all tie point labels
    for i=1:m
        nt=handles.ties(i,1);
        hobj1=findobj('String',nt,'Parent',handles.axes1,'Color','k');
        hobj2=findobj('String',nt,'Parent',handles.axes2,'Color','k');
        x=get(hobj1,'Position');
        delete(hobj1);
        delete(findobj('XData',x(1),'Parent',handles.axes1));
        y=get(hobj2,'Position');
        delete(findobj('XData',y(1),'Parent',handles.axes2));
        delete(hobj2);
    end
    % Sort tie points by position in signal
    if(~isempty(handles.ties))
        % Sort by core number
        [b,ind]=sort(handles.ties(:,2));
        j=-1;
        si=[];
        % Sort tie points within each core
        for i=0:max(handles.ties(:,2))
            subind=find(handles.ties(:,2)==i);
            if(~isempty(subind) & i~=j)
                [d,temp]=sort(handles.ties(subind,3));
                %size_temp=size(temp)
                si(end+1:end+size(subind))=subind(temp); %indices of sorted points
            end
        end
        %size_si=size(si)
        % Save sorted points with new label numbers
        handles.ties=[[1:m]' handles.ties(si,2:5)];
        guidata(gcbo,handles);
        plot_ties(h, eventdata, handles, varargin);
        
        % Test for crossed tie points
        j=-1;
        si=[];
        % Sort by core number
        for i=0:max(handles.ties(:,4))
            subind=find(handles.ties(:,4)==i);
            % Sort within each core
            if(~isempty(subind) & i~=j)
                [d,temp]=sort(handles.ties(subind,5));
                si(end+1:end+size(subind))=subind(temp);
                j=i;
            end
        end
        % Test if sorting of target points match sorting of signal points
        test=(handles.ties(:,5)~=handles.ties(si,5));
        if(sum(test)>0)
            disp('Warning: Crossed tie points:')
            index_numbers=find(test==1)'
        end
    end
else            % If gaps are selected
    [m,n]=size(handles.gaps);
    handles.ngaps=m+1;
    % Remove all gap labels from both axes
    delete(findobj('Parent',handles.axes1,'Color','r'));
    delete(findobj('Parent',handles.axes2,'Color','r'));

    if(~isempty(handles.gaps))
        d=[];
        e=[];
        handles.gaps(:,:);
        % find all gaps in signal
        ind1=find(handles.gaps(:,2)==1);
        [m,y]=size(ind1);
        % sort signal gaps
        if(m>0)
            c=handles.gaps(ind1,:);
            [b,ind]=sort(c(:,4));
            d=[[1:m]' c(ind,2:5)];  %sorted gaps with new labels
        end
        % find all gaps in target
        ind2=find(handles.gaps(:,2)==2);
        [n,y]=size(ind2);
        % sort target gaps
        if(n>0)
            c=handles.gaps(ind2,:);
            [b,ind]=sort(c(:,4));
            e=[[m+1:m+n]' c(ind,2:5)];   %sorted gaps with new labels
        end 
        % store modifications
        handles.gaps=[d;e];
        guidata(gcbo,handles);
        % replot gaps with new labels
        plot_ties(h, eventdata, handles, 1);
    end
end

% --------------------------------------------------------------------
function varargout = push_del_all_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_del_all.
% --------------------------------------------------------------------
% Delete all tie points or gaps from signal and target

handles=guidata(handles.figure1);
test=0;
if(nargin == 4)    % function called from other function (not button press)
    test=varargin{1};
end
% If tie points selected or indicated in argument list
if((get(handles.radio_tie,'Value')==1 & test==0) | test==1)
    % Delete all black objects (labels and markers from graphs
    ans=confirm(2);
    if(ans==1)
        delete(findobj('Parent',handles.axes1,'Color','k'));
        delete(findobj('Parent',handles.axes2,'Color','k'));
        % Clear tie points from memory
        handles.ties=[];
        handles.nties=1;
    end
end
% If gaps selected or indicated in argument list
if((get(handles.radio_tie,'Value')==0 & test==0) | test==2)
    ans=confirm(3);
    if(ans==1)
        % Delete all black objects (labels and markers from graphs
        delete(findobj('Parent',handles.axes1,'Color','r'));
        delete(findobj('Parent',handles.axes2,'Color','r'));
        % Clear tie points from memory
        handles.gaps=[];
        handles.ngaps=1;
    end
end
guidata(gcbo,handles);

% --------------------------------------------------------------------
function varargout = plot_ties(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.toggle_tie.
% --------------------------------------------------------------------
% Plot all tie points or gaps (called from other functions)
% Any varargin= plot gaps, empty varagin= plot tie points

% Get handles data
handles=guidata(handles.figure1);
% Read in signal and target
f1=get(handles.edit_sig,'String');
f2=get(handles.edit_targ,'String');
try
    s1=load(f1);
    [m n]=size(s1);
    xc1=n-1;
    yc1=n;
end
try
    s2=load(f2);
    [m n]=size(s2);
    xc2=n-1;
    yc2=n;
end

if(~isempty(varargin{1})) % Gaps
    gaps=handles.gaps;
    for i=1:size(gaps)
        if(gaps(i,2)==1)
            axes(handles.axes1)
        else
            axes(handles.axes2)
            s1=s2;
            xc1=xc2;
            yc1=yc2;
        end

        % Determine type of gap
        left=handles.gaps(i,4);
        right=handles.gaps(i,5);
        if(yc1==2)
            lend=min(s1(:,xc1));
            rend=max(s1(:,xc1));
        else
            core=handles.gaps(i,3);
            lend=min(s1(find(s1(:,1)==core),xc1));
            rend=max(s1(find(s1(:,1)==core),xc1));
        end
        if(left==right)
            pos=4;
        elseif(left==lend)
            if(right==rend)
                pos=6;
            else
                pos=1;
            end
        elseif(right==rend)
            pos=2;
        else
            pos=3;
        end
        ans=pos;
        label=num2str(handles.gaps(i,1));
        x1=left;
        x2=right;
        % Get data point positions for plots
        if(yc1==2)
            val=find(s1(:,xc1)<=x1);
            val2=find(s1(:,xc1)<=x2);
        else
            val=find(s1(:,1)==core & s1(:,xc1)<=x1);
            val2=find(s1(:,1)==core & s1(:,xc1)<=x2);
        end
       
        % Plot gaps and labels (squares indicate moveable ends)
        if(ans~=1)
            plot(x1,s1(val(end),yc1),'rs','MarkerSize',4);
            text(x1,s1(val(end),yc1),label,'VerticalAlignment','bottom',...
                'HorizontalAlignment','center','Color','r','FontSize',12);
        end
        if(ans==1 | ans==3 | ans==7)
            plot(x2,s1(val2(end),yc1),'rs','MarkerSize',4);
            text(x2,s1(val2(end),yc1),label,'VerticalAlignment','bottom',...
                'HorizontalAlignment','center','Color','r','FontSize',12);
        end
        % Signal/target plotted in red over gap area
        if(ans~=4)
            plot(s1(val(end):val2(end),xc1),s1(val(end):val2(end),yc1),'r');
        end
    end
            
else                     % Tie points indicated by empty varargin
    ties=handles.ties;
    for i=1:handles.nties-1
        % Signal side of tie point
        axes(handles.axes1)
        tie=ties(i,:);       % Info of selected tie point
        % Find data points for plotting
        if(xc1==1)
            val=find(s1(:,xc1)<=tie(3));
        else
            val=find(s1(:,xc1)<=tie(3) & s1(:,1)==tie(2));
            if(isempty(val))
                val=min(find(s1(:,xc1)>=tie(3) & s1(:,1)==tie(2)));
            end
        end            
        % Plot signal tie point and label
        plot(tie(3),s1(val(end),yc1),'k+','MarkerSize',5);
        text(tie(3),s1(val(end),yc1),num2str(tie(1)),'VerticalAlignment','bottom',...
            'HorizontalAlignment','center','FontSize',12);
        % Signal side of tie point
        axes(handles.axes2)
        if(xc2==1)
            val=find(s2(:,xc2)<=tie(5));
        else
            val=find(s2(:,xc2)<=tie(5) & s2(:,1)==tie(4));
        end          
        % Plot target tie point and label
        plot(tie(5),s2(val(end),yc2),'k+','MarkerSize',5);
        text(tie(5),s2(val(end),yc2),num2str(tie(1)),'VerticalAlignment','bottom',...
            'HorizontalAlignment','center','FontSize',12);
    end
end

% --------------------------------------------------------------------
function varargout = end_pen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.end_pen.
% Text box for "Endpoint" (aka "No Match") penalty

% --------------------------------------------------------------------
function varargout = speed_pen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.speed_pen.
% Text box for "Speed" penalty (penalizes speeds away from target speed)

% --------------------------------------------------------------------
function varargout = spch_pen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.spch_pen.
% Text box for "Speed Change" penalty (penalizes changes in matching speed)

% --------------------------------------------------------------------
function varargout = tie_pen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.tie_pen.
% Text box for "Tie Point" penalty (penalizes non-matching tie points)

% --------------------------------------------------------------------
function varargout = gap_pen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.gap_pen.
% Text box for "Gap" penalty (penalizes matched gap sizes not equal original size)

% --------------------------------------------------------------------
function varargout = toggle_auto_pen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.toggle_auto_pen.
% --------------------------------------------------------------------
% Tries to calculate reasonable penalty values based on magnitude and
% similarity of signal and target (Mostly an order-of-magnitude guide,
% and may be significantly off in some situations).

if(get(handles.toggle_auto_pen,'Value')==1)
    % Load signal and target 
    f1=get(handles.edit_sig,'String');
    f2=get(handles.edit_targ,'String');
    try
        s1=load(f1);
        [m n]=size(s1);
        xc1=n-1;
        yc1=n;
    end
    try
        s2=load(f2);
        [m n]=size(s2);
        xc2=n-1;
        yc2=n;
    end
    
    % Find subset of signal & target included in match
    st1=str2num(get(handles.edit_beg1,'String'));
    e1=str2num(get(handles.edit_end1,'String'));
    sub1=s1(find(s1(:,xc1)>=st1 & s1(:,xc1)<=e1),:);
    st2=str2num(get(handles.edit_beg2,'String'));
    e2=str2num(get(handles.edit_end2,'String'));
    sub2=s2(find(s2(:,xc2)>=st2 & s2(:,xc2)<=e2),:);
    
    % Calculate means and standard deviations
    mean1=mean(sub1(:,yc1));
    mean2=mean(sub2(:,yc2));
    m=2*abs(mean1-mean2);
    std1=std(sub1(:,yc1));
    std2=std(sub2(:,yc2));
    s=std1;
    if(std2 > s)
        s=std2;
    end
    s=s^2;       % largest variance from signal or target
    d=e2-st2;    % length of target included in match
    
    set(handles.end_pen,'String',num2str(round(150*(s+m^2))/10));
    set(handles.speed_pen,'String',num2str(round(350*(.2*s+m^2))/100));
    set(handles.spch_pen,'String',num2str(round(350*(.15*s+m^2))/100));
    set(handles.tie_pen,'String',num2str(round(5000*s/d+m^2)));
    set(handles.gap_pen,'String',num2str(round(1000*s/d+.8*m^2)));
end

% --------------------------------------------------------------------
function answer = speeds2str(handles)
% --------------------------------------------------------------------
% Translates speed arrays into strings for configuration file or display
% and updates prefered speed pop-up menu

[m n]=size(handles.speed_num);             % Find number of speeds
t=sscanf(handles.target_speed,'%d : %d');  % Prefered match speed

s='';
for i=1:n
    s(end+1)=num2str(handles.speed_num(i));   % Numerator of speed ratio 
    s(end+1)=':';
    s(end+1)=num2str(handles.speed_den(i));   % Denominator of speed ratio 
    if(t(1)==handles.speed_num(i) & t(2)==handles.speed_den(i))
        % Keep prefered speed display indexed correctly
        set(handles.popupmenu1,'Value',i);    
    end
    if(i~=n)
        s(end+1)=',';
    end
end
% If prefered speed is no longer present, set new prefered speed to first speed in list
if(n<get(handles.popupmenu1,'Value'))
    set(handles.popupmenu1,'Value',1);
end
% Set speeds to display in pop-up menu for user to select prefered speed
pops=s;
pops(find(pops==','))='|';
set(handles.popupmenu1,'String',pops);
guidata(handles.figure1,handles)
% Return string of speeds 
answer=s;

% --------------------------------------------------------------------
function varargout = button_get_speeds_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.button_get_speeds.
% --------------------------------------------------------------------
% Button for user to set new match speeds: calls getspeed.m; translates string
% into array and calls speeds2str to update prefered speed pop-up menu

speeds=getspeed(handles.speed_num,handles.speed_den);
ind=find(speeds(:)==':');
ind2=find(speeds(:)==',' | speeds(:)==';' | speeds(:)=='[' | speeds(:)=='(');
speeds(ind2)=' ';
[m n]=size(ind);
num=[];
den=[];
num(1)=sscanf(speeds(1:ind(1)-1),'%d');
for i=1:m-1
    b=sscanf(speeds(ind(i)+1:ind(i+1)-1),'%d %d');
    den(end+1)=b(1);
    num(end+1)=b(2);
end
den(end+1)=sscanf(speeds(ind(m)+1:end),'%d');

handles.speed_num=num;
handles.speed_den=den;
guidata(gcbo,handles);

speeds=speeds2str(handles);

% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu1.
% --------------------------------------------------------------------
% Pop-up menu for user to select prefered ("target") matching speed from
% the list of possible speeds

val=get(handles.popupmenu1,'Value');
sp='';
sp(1)=num2str(handles.speed_num(val));
sp(2)=':';
sp(3)=num2str(handles.speed_den(val));
handles.target_speed=sp;
guidata(gcbo,handles);

% --------------------------------------------------------------------
function varargout = edit_axeslim1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_axeslim1.
% --------------------------------------------------------------------
% Text box to set signal plot limits  
% Text format: x-min, x-max, y-min, y-max separated by commas, semi-colons, or spaces

axes(handles.axes1)
lim1=get(handles.edit_axeslim1,'String');
ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
lim1(ind)=' ';
b=sscanf(lim1,'%g %g %g %g');
if(max(size(b))==4 & b(2)>b(1) & b(4)>b(3))
    xlim(handles.axes1, [b(1) b(2)]);
    ylim(handles.axes1, [b(3) b(4)]);
    if(b(1)>get(handles.slider1,'Max'))
        set(handles.slider1,'Value',get(handles.slider1,'Max'));
    elseif(b(1)<get(handles.slider1,'Min'))    
        set(handles.slider1,'Value',get(handles.slider1,'Min'));
    else
        set(handles.slider1,'Value',b(1),'SliderStep',[.03 .06]);
    end
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];

    if(get(handles.transfer_axeslim,'Value')==1)
        xlim(handles.axes2, [b(1) b(2)]);
        ylim(handles.axes2, [b(3) b(4)]);
        if(b(1)>get(handles.slider2,'Max'))
            set(handles.slider2,'Value',get(handles.slider2,'Max'));
        elseif(b(1)<get(handles.slider2,'Min'))    
            set(handles.slider2,'Value',get(handles.slider2,'Min'));
        else
            set(handles.slider2,'Value',b(1),'SliderStep',[.03 .06]);
        end
        set(handles.edit_axeslim2,'String',c);
        if(handles.match_fig~=0)
            h=handles.results_axes;
            for i=1:max(size(h))
                xlim(handles.results_axes(i), [b(1) b(2)]);
            end
            ylim(handles.results_axes(1), [b(3) b(4)]);
            set(handles.results_lim,'String',c);
        end
    end
else
    b=axis;
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim1,'String',c);
end


% --------------------------------------------------------------------
function varargout = edit_axeslim2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit_axeslim2.
% --------------------------------------------------------------------
% Text box to set target plot limits  
% Text format: x-min, x-max, y-min, y-max separated by commas, semi-colons, or spaces

axes(handles.axes2)
lim1=get(handles.edit_axeslim2,'String');
ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
lim1(ind)=' ';
b=sscanf(lim1,'%g %g %g %g');
if(max(size(b))==4 & b(2)>b(1) & b(4)>b(3))
    xlim(handles.axes2, [b(1) b(2)]);
    ylim(handles.axes2, [b(3) b(4)]);
    if(b(1)>get(handles.slider2,'Max'))
        set(handles.slider2,'Value',get(handles.slider2,'Max'));
    elseif(b(1)<get(handles.slider2,'Min'))    
        set(handles.slider2,'Value',get(handles.slider2,'Min'));
    else
        set(handles.slider2,'Value',b(1),'SliderStep',[.03 .06]);
    end
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];

    if(get(handles.transfer_axeslim,'Value')==1)
        xlim(handles.axes1, [b(1) b(2)]);
        ylim(handles.axes1, [b(3) b(4)]);
        if(b(1)>get(handles.slider1,'Max'))
            set(handles.slider1,'Value',get(handles.slider1,'Max'));
        elseif(b(1)<get(handles.slider1,'Min'))    
            set(handles.slider1,'Value',get(handles.slider1,'Min'));
        else
            set(handles.slider1,'Value',b(1),'SliderStep',[.03 .06]);
        end
        set(handles.edit_axeslim1,'String',c);
        if(handles.match_fig~=0)
            h=handles.results_axes;
            for i=1:max(size(h))
                xlim(handles.results_axes(i), [b(1) b(2)]);
            end
            ylim(handles.results_axes(1), [b(3) b(4)]);
            set(handles.results_lim,'String',c);
        end
    end
else
    b=axis;
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim2,'String',c);
end



% --------------------------------------------------------------------
function varargout = results_lim_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.results_lim.
% --------------------------------------------------------------------
% Text box to set results plot limits  
% Text format: x-min, x-max, y-min, y-max separated by commas, semi-colons, or spaces

if(ishandle(handles.match_fig))
    lim1=get(handles.results_lim,'String');
    ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
    lim1(ind)=' ';
    b=sscanf(lim1,'%g %g %g %g');
    h=handles.results_axes;
    if(max(size(b))==4 & b(2)>b(1) & b(4)>b(3))
        for i=1:max(size(h))
            xlim(handles.results_axes(i), [b(1) b(2)]);
        end
        ylim(handles.results_axes(1), [b(3) b(4)]);
        
        if(get(handles.transfer_axeslim,'Value')==1)
            c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
            xlim(handles.axes1, [b(1) b(2)]);
            ylim(handles.axes1, [b(3) b(4)]);
            if(b(1)>get(handles.slider1,'Max'))
                set(handles.slider1,'Value',get(handles.slider1,'Max'));
            elseif(b(1)<get(handles.slider1,'Min'))    
                set(handles.slider1,'Value',get(handles.slider1,'Min'));
            else
                set(handles.slider1,'Value',b(1),'SliderStep',[.03 .06]);
            end
            set(handles.edit_axeslim1,'String',c);
            xlim(handles.axes2, [b(1) b(2)]);
            ylim(handles.axes2, [b(3) b(4)]);
            if(b(1)>get(handles.slider2,'Max'))
                set(handles.slider2,'Value',get(handles.slider2,'Max'));
            elseif(b(1)<get(handles.slider2,'Min'))    
                set(handles.slider2,'Value',get(handles.slider2,'Min'));
            else
                set(handles.slider2,'Value',b(1),'SliderStep',[.03 .06]);
            end
            set(handles.edit_axeslim2,'String',c);
        end            
    else
        b=axis(handles.results_axes(1));
        c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
        set(handles.results_lim,'String',c);
    end
end
  
% --------------------------------------------------------------------
function varargout = transfer_axeslim_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.transfer_axeslim.
% --------------------------------------------------------------------
% Set signal and results axes to match target axes (Other functions will maintain 
% synchronization while button is depressed.)

if(get(handles.transfer_axeslim,'Value')==1)
    lim1=get(handles.edit_axeslim2,'String');
    ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
    lim1(ind)=' ';
    b=sscanf(lim1,'%g %g %g %g');
    xlim(handles.axes1, [b(1) b(2)]);
    ylim(handles.axes1, [b(3) b(4)]);
    if(b(1)>get(handles.slider1,'Max'))
        set(handles.slider1,'Value',get(handles.slider1,'Max'));
    elseif(b(1)<get(handles.slider1,'Min'))    
        set(handles.slider1,'Value',get(handles.slider1,'Min'));
    else
        set(handles.slider1,'Value',b(1),'SliderStep',[.03 .06]);
    end
    c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
    set(handles.edit_axeslim1,'String',c);

    if(handles.match_fig~=0)
        h=handles.results_axes;
        for i=1:max(size(h))
            xlim(handles.results_axes(i), [b(1) b(2)]);
        end
        ylim(handles.results_axes(1), [b(3) b(4)]);
        set(handles.results_lim,'String',c);
    end
end


% --------------------------------------------------------------------
function varargout = flipy_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.flipy.
% --------------------------------------------------------------------
% When depressed y-axis of signal is reversed

if(get(handles.flipy,'Value')==1)
    set(handles.axes1,'YDir','reverse');
else
    set(handles.axes1,'YDir','normal');
end

% --------------------------------------------------------------------
function varargout = flipx_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.flipx.
% --------------------------------------------------------------------
% When depressed x-axis of signal is reversed
if(get(handles.flipx,'Value')==1)
    set(handles.axes1,'XDir','reverse');
else
    set(handles.axes1,'XDir','normal');
end

% --------------------------------------------------------------------
function varargout = flipy2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.flipy.
% --------------------------------------------------------------------
% When depressed y-axis of target is reversed
if(get(handles.flipy2,'Value')==1)
    set(handles.axes2,'YDir','reverse');
else
    set(handles.axes2,'YDir','normal');
end

% --------------------------------------------------------------------
function varargout = flipx2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.flipx.
% --------------------------------------------------------------------
% When depressed x-axis of target is reversed
if(get(handles.flipx2,'Value')==1)
    set(handles.axes2,'XDir','reverse');
else
    set(handles.axes2,'XDir','normal');
end


% --------------------------------------------------------------------
function varargout = find_sig_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
% --------------------------------------------------------------------
% User scans through a directory list to select the file containing the signal
sig=getconf(6);
if(~strcmp(sig,'cancel'))
    set(handles.edit_sig,'String',sig);
    edit_sig_Callback(h, eventdata, handles);
end
    
% --------------------------------------------------------------------
function varargout = find_targ_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.find_sig.
% --------------------------------------------------------------------
% User scans through a directory list to select the file containing the target
targ=getconf(6);
if(~strcmp(targ,'cancel'))
    set(handles.edit_targ,'String',targ);
    edit_targ_Callback(h, eventdata, handles);
end


% --------------------------------------------------------------------
function varargout = load_match_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.load_mactch.
% --------------------------------------------------------------------
% Plots the signal aligned to the target and its relative sedimentation rate
% in a new figure window using the result files generated by Match
% (The configuration file must be generated and run in Match before this will work.)
% Will plot all signal alignments in multi-signal mode.

%Open new window
match_fig=figure;
handles.match_fig=match_fig;
guidata(gcbo,handles);
%Set window's position on screen
pos=get(match_fig,'Position');
pos(1)=pos(1)-pos(3)/5;
pos(3)=pos(3)*1.4;
set(match_fig,'Position',pos);

% Read in match file
name=[get(handles.edit_out,'String') '.match'];
set(match_fig,'Name',name);
xmatch=load(name,'-ascii');

% Get signal filename(s)
if(get(handles.multi_sig,'Value')==0)
    num=1;
else
    handles=guidata(handles.figure1);
    sigs=handles.moresig;
    [m,n]=size(sigs);
    num=m/2+1;
end

% Loop through number of signals
for i=1:num
    % Load target and aligned signal
    if(i==1)
        sig=get(handles.edit_sig,'String');
        targ=get(handles.edit_targ,'String');
    else
        sig=deblank(sigs(2*(i-1)-1,:));
        targ=deblank(sigs(2*(i-1),:));
    end
    f1=[sig '.new'];
    f2=[targ];
    n1=0;
    n2=0;
    m1=[];
    m2=[];
    try
        s1=load(f1);
        [m n1]=size(s1);
    catch
        disp('Error: Invalid file')
        f1
    end
    try
        s2=load(f2);
        [m n2]=size(s2);
    catch
        disp('Error: Invalid file')
        f2
    end
    
    if(n1~=0 & n2~=0)
        h(i)=subplot(num+1,1,i);
        hold on;
        % Find gaps in data
        [m,p]=find(isnan(s1));
        j1=size(m);
        if(n1==3 & isempty(m))
            [m,p]=find(s1(1:end-1,1)~=s1(2:end,1));
            j1=size(m);
            for k=1:j1
                s1(m(k)+k:end+1,:)=[NaN NaN NaN;s1(m(k)+k:end,:)];
                m(k)=m(k)+k-1;
            end
            m1=m;
        elseif(~isempty(m))
            m1=m-1;
        end
        
        [m,p]=find(isnan(s2));
        j2=size(m);
        if(n2==3 & isempty(m))
            [m,p]=find(s2(1:end-1,1)~=s2(2:end,1));
            j2=size(m);
            for k=1:j2
                s2(m(k)+k:end+1,:)=[NaN NaN NaN;s2(m(k)+k:end,:)];
                m(k)=m(k)+k-1;
            end
            m2=m;
        elseif(~isempty(m))
            m2=m-1;
        end

        % Plot aligned signal and target
        plot(s2(:,n2-1),s2(:,n2),'Color',[.35 .35 1]);
        plot(s1(:,n1-1),s1(:,n1),'m');
        plot(s1(m1(:),n1-1),s1(m1(:),n1),'.','Color',[1 .01 .01],'MarkerSize',9);
        plot(s2(m2(:),n2-1),s2(m2(:),n2),'.','Color',[1 .01 .01],'MarkerSize',9);
        t=[f1(1:end-4) ' aligned to ' f2];
        title(t)
        %legend(f2,f1,0);
        
        % Set axes limits
        if(get(handles.transfer_axeslim,'Value')==1)
            lim1=get(handles.edit_axeslim2,'String');
            ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
            lim1(ind)=' ';
            b=sscanf(lim1,'%g %g %g %g');
            c=[num2str(b(1)),', ',num2str(b(2)),', ',num2str(b(3)),', ',num2str(b(4))];
            set(handles.results_lim,'String',c);
            xlim([b(1) b(2)]);
            if(i==1)
                ylim([b(3) b(4)]);
            end
        else
            xlim([min(s1(:,n1-1)) max(s1(:,n1-1))]);
        end
        xl=xlim;
        % Set axes orientation to be the same as the target axes
        if(get(handles.flipy2,'Value')==1)
            set(gca,'YDir','reverse');
        end
        if(get(handles.flipx2,'Value')==1)
            set(gca,'XDir','reverse');
        end

        % Plot tie points
        ties=handles.ties;
        [m,n]=size(ties);
        for i=1:m
            tie=ties(i,:);
            if(n1==2)
                newx=interp1(xmatch(:,2),xmatch(:,4),tie(3));
                val=find(s1(:,n1-1)<=newx);
                y=s1(val(end),n1);
                plot(newx,y,'k+','MarkerSize',5);
                text(newx,y,num2str(tie(1)),'Color','k','VerticalAlignment','bottom',...
                    'HorizontalAlignment','center','FontSize',12);
            else
                try
                    ind=find(s1(:,1)==tie(2));
                    if(~isempty(ind))
                        newx=interp1(xmatch(ind,2),xmatch(ind,4),tie(3));
                        val=find(s1(ind,2)<=newx);
                        y=s1(val(end),n1);
                        plot(newx,y,'k+','MarkerSize',5);
                        text(newx,y,num2str(tie(1)),'Color','k','VerticalAlignment','bottom',...
                            'HorizontalAlignment','center','FontSize',12);
                    end
                end
            end            
            if(n2==2)
                val=find(s2(:,n2-1)<=tie(5));
            else
                val=find(s2(:,n2-1)<=tie(5));  % & s1(:,1)==tie(4));
            end            
            y=s2(val(end),n2);
            plot(tie(5),y,'k+','MarkerSize',5);
            text(tie(5),y,num2str(tie(1)),'Color','k','VerticalAlignment','top',...
                'HorizontalAlignment','center','FontSize',12);
        end;
        % Plot gaps
        gaps=handles.gaps;
        [m,n]=size(gaps);
        for i=1:m
            if(gaps(i,2)==1)
                try
                    newx=interp1(xmatch(:,2),xmatch(:,4),gaps(i,3));
                    val=find(s1(:,n1-1)<=newx);
                    y=s1(val(end),n1);
                    plot(newx,y,'rs','MarkerSize',4);
                    text(newx,y,num2str(gaps(i,1)),'Color','r','VerticalAlignment','bottom',...
                        'HorizontalAlignment','center','FontSize',12);
                end
            else
                val=find(s2(:,n2-1)<=gaps(i,3));
                y=s2(val(end),n2);
                plot(gaps(i,3),y,'rs','MarkerSize',4);
                text(gaps(i,3),y,num2str(gaps(i,1)),'Color','r','VerticalAlignment','top',...
                    'HorizontalAlignment','center','FontSize',12);
            end    
        end
    end    
end    

% Plot relative accumulation (sedimentation) rate 
t=xmatch(:,4);
d=xmatch(:,2);
srate=d(1:(size(t,1)-1));
for i=1:(size(t,1)-1)
    srate(i)=(d(i+1)-d(i))/(t(i+1)-t(i));
end
h(num+1)=subplot(num+1,1,num+1);
plot(t(1:(size(t,1)-1)),srate,'b-');
xlim(xl)
if(get(handles.flipx2,'Value')==1)
    set(gca,'XDir','reverse');
end
title('Relative accumulation rate')

handles.results_axes=h;
guidata(handles.figure1,handles);    


% --------------------------------------------------------------------
function varargout = multi_sig_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.multi_sig.
% --------------------------------------------------------------------
% Allows user to enter up to three additional signals and targets to be
% matched simultaneously. (All signals must have same depth/time range, and
% all targets must have same depth/time range.)

if(get(handles.multi_sig,'Value')==1)
    files=get_multi_sig;
    if(isempty(files))
        set(handles.multi_sig,'Value',0);
    else
        j=1;
        for i=1:6
            q=isspace(files(i,:));
            if(min(q)==0)
                f(j,:)=char(files(i,:));
                j=j+1;
            end
        end
        handles.moresig=f;
        guidata(gcbo,handles);
        plot_moresig(h, eventdata, handles);
    end
end


% --------------------------------------------------------------------
function varargout = refresh_moresig_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.refresh_moresig.
% --------------------------------------------------------------------
% Replots the multi-signal plot to reflect any changes in the axes limits
if(get(handles.multi_sig,'Value')==1)
    plot_moresig(h, eventdata, handles);
end

% --------------------------------------------------------------------
function varargout = plot_moresig(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
% Opens an additional window that plots all signals and targets for 
% multi-signal matching. (All axes, tie point, and gap changes must be 
% performed on original signal and target axes. Pressing the "Replot"
% button will update the multi-signal window to reflect these changes.

handles=guidata(handles.figure1);
sigs=handles.moresig;
[m,n]=size(sigs);
num=m/2+1;
width=.9/(m+2);

if(handles.sig_fig==0)
    sf=0;
    sig_fig=figure;
    pos=get(sig_fig,'Position');
    pos(2)=.3*pos(2);
    pos(4)=pos(4)*1.5;
    set(sig_fig,'Position',pos);
	set(sig_fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    set(sig_fig,'Name','All Signals');
    handles.sig_fig=sig_fig;
    guidata(handles.figure1,handles);
else
    sf=1;
    sig_fig=handles.sig_fig;
    figure(sig_fig)
    clf
end

for i=1:num
    if(i==1)
        sig=get(handles.edit_sig,'String');
        targ=get(handles.edit_targ,'String');
    else
        sig=deblank(sigs(2*(i-1)-1,:));
        targ=deblank(sigs(2*(i-1),:));
    end
    yc1=0;
    yc2=0;
    try
        s1=load(sig);
        [m yc1]=size(s1);
    catch
        disp('Error: Invalid signal file')
        sig
    end
    try
        s2=load(targ);
        [m yc2]=size(s2);
    catch
        disp('Error: Invalid target file')
        targ
    end
    
    if(yc1~=0 & yc2~=0)
        s1=load(sig);
        s2=load(targ);
        subplot('Position',[.05+width*2*(i-1) .08 width-.005 .84]);

        [m,p]=find(isnan(s1));
        j=size(m);
        if(yc1==3 & isempty(m))
            [m,p]=find(s1(1:end-1,1)~=s1(2:end,1));
            j=size(m);
            for k=1:j
                s1(m(k)+k:end+1,:)=[NaN NaN NaN;s1(m(k)+k:end,:)];
                m(k)=m(k)+k-1;
            end
        elseif(~isempty(m))
            m=m-1;
        end
        
        plot(s1(:,yc1),s1(:,yc1-1),'m');
        hold on;
        for k=1:j
            plot(s1(m(k),yc1),s1(m(k),yc1-1),'.','Color',[1 .01 .01],'MarkerSize',9);
        end
        axis ij;
        gca1=gca;
        set(gca,'XTickLabelMode','manual');
        if(i==1)
            st1=str2num(get(handles.edit_beg1,'String'));
            e1=str2num(get(handles.edit_end1,'String'));
            ylim1=[st1 e1];
        else
            set(gca,'YTickLabelMode','manual');
        end
        if(sf~=0)
            ylim1=xlim(handles.axes1);
            set(gca,'YLim',ylim1);
        end
        title(sig)
        if(get(handles.transfer_axeslim,'Value')==1)
            lim1=get(handles.edit_axeslim2,'String');
            ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
            lim1(ind)=' ';
            b=sscanf(lim1,'%g %g %g %g');
            ylim([b(1) b(2)]);
            if(i==1)
                xlim([b(3) b(4)]);
            end
        end
        hold on;

        subplot('Position',[.045+width*(2*(i-1)+1) .08 width-.005 .84]);
        [m,p]=find(isnan(s2));
        j=size(m);
        if(yc2==3 & isempty(m))
            [m,p]=find(s2(1:end-1,1)~=s2(2:end,1));
            j=size(m);
            for k=1:j
                s2(m(k)+k:end+1,:)=[NaN NaN NaN;s2(m(k)+k:end,:)];
                m(k)=m(k)+k-1;
            end
        elseif(~isempty(m))
            m=m-1;
        end
           
        set(gca,'YTickLabelMode','manual');
        plot(s2(:,yc2),s2(:,yc2-1),'Color',[.35 .35 1]);
        hold on;
        for k=1:j
            plot(s2(m(k),yc2),s2(m(k),yc2-1),'.','Color',[1 .01 .01],'MarkerSize',9);
        end
        set(gca,'XTickLabelMode','manual');
        set(gca,'YTickLabelMode','manual');

        axis ij;
        gca2=gca;
        if(i==1)
            st2=str2num(get(handles.edit_beg2,'String'));
            e2=str2num(get(handles.edit_end2,'String'));
            ylim2=[st2 e2];
        end
        if(sf~=0)
            ylim2=xlim(handles.axes2);
            set(gca,'YLim',ylim2);
        end
        if(i==num)
            set(gca,'YAxisLocation','right');
            set(gca,'YTickLabelMode','auto');
        else
            set(gca,'YTickLabelMode','manual');
        end
        title(targ)
        if(get(handles.transfer_axeslim,'Value')==1)
            lim1=get(handles.edit_axeslim2,'String');
            ind=find(lim1(:)==',' | lim1(:)==';' | lim1(:)=='[' | lim1(:)=='(');
            lim1(ind)=' ';
            b=sscanf(lim1,'%g %g %g %g');
            ylim([b(1) b(2)]);
            if(i==1)
                xlim([b(3) b(4)]);
            end
        end
        hold on;
    end
    
    xc1=yc1-1;
    xc2=yc2-1;
    gaps=handles.gaps;
    
    for ii=1:size(gaps)
        if(gaps(ii,2)==1)
            axes(gca1)
            s=s1;
            xc=xc1;
            yc=yc1;
        else
            axes(gca2)
            s=s2;
            xc=xc2;
            yc=yc2;
        end

        left=handles.gaps(ii,4);
        right=handles.gaps(ii,5);
        if(yc==2)
            lend=min(s(:,xc));
            rend=max(s(:,xc));
        else
            core=handles.gaps(ii,3);
            lend=min(s(find(s(:,1)==core),xc));
            rend=max(s(find(s(:,1)==core),xc));
        end
        if(left==right)
            pos=4;
        elseif(left==lend)
            if(right==rend)
                pos=6;
            else
                pos=1;
            end
        elseif(right==rend)
            pos=2;
        else
            pos=3;
        end
        ans=pos;
        label=num2str(handles.gaps(ii,1));

        x1=left;
        x2=right;
        if(yc==2)
            val=find(s(:,xc)<=x1);
            val2=find(s(:,xc)<=x2);
        else
            val=find(s(:,1)==core & s(:,xc)<=x1);
            val2=find(s(:,1)==core & s(:,xc)<=x2);
        end
       
        if(ans~=1)
            plot(s(val(end),yc),x1,'rs','MarkerSize',4);
            text(s(val(end),yc),x1,label,'VerticalAlignment','middle',...
                'HorizontalAlignment','left','Color','r','FontSize',12);
        end
        if(ans==1 | ans==3 | ans==7)
            plot(s(val2(end),yc),x2,'rs','MarkerSize',4);
            text(s(val2(end),yc),x2,label,'VerticalAlignment','middle',...
                'HorizontalAlignment','left','Color','r','FontSize',12);
        end
        if(ans~=4)
            plot(s(val(end):val2(end),yc),s(val(end):val2(end),xc),'r');
        end
    end
            
    ties=handles.ties;
    for ii=1:handles.nties-1
        axes(gca1)
        tie=ties(ii,:);
        if(xc1==1)
            val=find(s1(:,xc1)<=tie(3));
        else
            val=find(s1(:,xc1)<=tie(3) & s1(:,1)==tie(2));
        end    
        if(~isempty(val))
            plot(s1(val(end),yc1),tie(3),'k+','MarkerSize',5);
            text(s1(val(end),yc1),tie(3),num2str(tie(1)),'VerticalAlignment','middle',...
                'HorizontalAlignment','right','FontSize',12);
        end
        axes(gca2)
        if(xc2==1)
            val=find(s2(:,xc2)<=tie(5));
        else
            val=find(s2(:,xc2)<=tie(5) & s2(:,1)==tie(4));
        end            
        if(~isempty(val))
            plot(s2(val(end),yc2),tie(5),'k+','MarkerSize',5);
            text(s2(val(end),yc2),tie(5),num2str(tie(1)),'VerticalAlignment','middle',...
                'HorizontalAlignment','right','FontSize',12);
        end 
    end
end 

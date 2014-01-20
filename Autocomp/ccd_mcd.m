function comphandles = ccd_mcd(varargin)
% CCD_MCD Application M-file for ccd_mcd.fig
%    FIG = CCD_MCD launch ccd_mcd GUI.
%    CCD_MCD('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 08-Jun-2006 10:56:19

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

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
	guidata(fig, handles);
    axes(handles.axes3)
    %a=handles.axes3;
    
    H=varargin{1};
    xlabel('MCD-MBSF offset (m)     ')
    set(handles.radio_ccdaxis,'String','MCD');
    set(handles.radio_mcdaxis,'String','CMCD');
    set(handles.radio_mcdaxis,'Position',[.259 .015 .116 .025]);
    set(handles.save_filename,'String',get(H.cdsname,'String'));
    H.apply=0;
    H.test=1;
    mcdtie=H.mcdtie;
    ind=find(mcdtie(:,1)==0);
    if isempty(ind)
        handles.newtie=[];
    else
        handles.newtie=mcdtie(ind,:);
    end        
    ind=find(mcdtie(:,1)~=0);
    if isempty(ind)
        handles.mcdtie_no_offset=[];
    else
        handles.mcdtie_no_offset=mcdtie(ind,:);
    end        
    handles.topoff=0;
    handles.H=H;
    handles.copyH=H;
	guidata(fig, handles);
    
    comp=H.comp;
    set(handles.targ_start,'String',num2str(comp(1,5)));
    set(handles.targ_end,'String',num2str(comp(end,6)));
    if(~isempty(H.eqn))
        polyorder=num2str(max(size(H.eqn))-1);
        s=['Polynomial fit has order ' polyorder];
        disp(s)
    else
        set(handles.radio_linear,'Value',1);
        set(handles.radio_poly,'Value',0);
        disp('Offset is linearly interpolated from mcd tie points')
    end
    if(~isempty(H.top))
        set(handles.mcd_start,'String',num2str(H.top));
        set(handles.targ_start,'String',num2str(H.start));
        set(handles.top_fixed,'Value',1);
        set(handles.top_auto,'Value',0);
    end
    figure(fig)
    plot_offset(fig, handles)
    
    uiwait(fig);
    %waitfor(fig,'BeingDeleted','on')
    if ~ishandle(fig)
        comphandles=H;
    else
        handles=guidata(fig);
        if(handles.apply==1)
            H=handles.H;
            H.apply=handles.apply;
            H.mcdtie=[handles.mcdtie_no_offset; handles.newtie];
            if(get(handles.top_fixed,'Value')==1)
                H.top=str2num(get(handles.mcd_start,'String'));
                H.start=str2num(get(handles.targ_start,'String'));
            else
                H.top=[];
                H.start=[];
            end
            comphandles = H;
        else
            comphandles=H;
        end
        delete(fig);
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
function varargout = plot_offset(h, handles, varargin)
% Stub for Callback of the uicontrol handles.toggle_tie.

%disp('Plot')
% if(~isempty(varargin))
%     varargin{1}
% end
handles=guidata(handles.mcdfigure);
figure(handles.mcdfigure);
a=handles.axes3;
H=handles.H;
axes(a);

%disp('plot_offset')

%mcdtie= hole, mbsf/mcd, ccd
offset=H.ccdoffset;
mcdtie=H.mcdtie;
num=H.numholes;

if(get(handles.radio_ccdaxis,'Value')==1)
    plot(offset(:,2),offset(:,1),'b')
    ylabel('MCD')
else
    plot(offset(:,2),offset(:,1)-offset(:,2),'k')
    ylabel('CMCD')
end
    hold on

for j=1:num
    if(round(j/2)==j/2)
        gr=.5;
    else
        gr=0;
    end
    color=[0 gr (j-1)/(num-1)];
    if(~isempty(mcdtie))
        ind=find(mcdtie(:,1)==j);
        if(get(handles.radio_ccdaxis,'Value')==1)
            plot(mcdtie(ind,3)-mcdtie(ind,2),mcdtie(ind,3),'+',...
                'Color',color,'LineWidth',2,'MarkerSize',5)
        else
            plot(mcdtie(ind,3)-mcdtie(ind,2),mcdtie(ind,2),'+',...
                'Color',color,'LineWidth',2,'MarkerSize',5)
        end
    end
end
if(~isempty(mcdtie))
    ind=find(mcdtie(:,1)==0);
    if(get(handles.radio_ccdaxis,'Value')==1)
        plot(mcdtie(ind,3)-mcdtie(ind,2),mcdtie(ind,3),'rx','LineWidth',2)
    else
        plot(mcdtie(ind,3)-mcdtie(ind,2),mcdtie(ind,2),'rx','LineWidth',2)
    end
end
xlabel('MCD-MBSF offset')
xlabel('MCD-MBSF offset (m)     ')
axis tight
mn=get(handles.targ_start,'String');
mn=str2num(mn);
mx=get(handles.targ_end,'String');
mx=str2num(mx);
if(~isempty(mx) & ~isempty(mn))
    ylim([mn mx])
end
axis ij
hold off

%---------------------------------------------------------------------------------------
function varargout = line_interp(h, handles, varargin)
% Linearly interpolate between tie points to convert from ccd to coretop mbsf

H=handles.H;
ccd=H.ccd;
comp=H.comp;
newtie=handles.newtie;
answer=confirm(7);
if(answer==1)   % Use coretops for hole 1 and user-defined tie points
    ct=handles.mcdtie_no_offset;
    %mcdtie= hole, mbsf/mcd, ccd
    ind=find(ct(:,1)==1); 
    if(~isempty(ind))
        oldct=ct(ind,:);
        ct=[ct(ind,:);newtie];
    else
        oldct=[];
        ct=newtie;
    end
else            % Use only user-defined mcd tie points
    ct=newtie;
    if(~isempty(ct))
        ind=find(ct(:,1)==0);
        newtie=newtie(ind,:);
    end
    oldct=[];
end

topoff=0;
top=ccd(1);
if(get(handles.top_fixed,'Value')==1)
    top=str2num(get(handles.mcd_start,'String'));
    if(~isempty(top))
        topoff=ccd(1)-top;
    end
end

if(isempty(ct))
    H.ccdoffset=[ccd(1) topoff; max(ccd) topoff];
    H.eqn=topoff;
    handles.H=H;
    guidata(h,handles);
    return
end

mx=max(ct(:,3));   % last mcdtie, excluding target core tops
if(isempty(mx))
    mx=max(ct(:,3));
end
%x1=[0:2:mx]';           % x values are in cumulative depth
ct=ct(find(ct(:,3)<=1.2*mx),:);
[b,i,j]=unique(ct(:,3));          % indices which exclude any possible duplicates
sz=max(size(i));

if(~isempty(i))
    ct=ct(i,:); 
    %mcdtie= hole, mbsf/mcd, ccd
    [y,ind]=sort(ct(:,3));
    [y,jnd]=sort(ct(:,2));
    if(~isempty(find(ind~=jnd)))
        disp('ERROR: Crossed mcd tie points. Cannot interpolate.')
        H.ccdoffset=[ccd(1) topoff; max(ccd) topoff];
        H.eqn=topoff;
        handles.H=H;
        guidata(h,handles);
        return
    end
end    

ct=ct(ind,:);
    
if(ccd(1)<ct(1,3))
    ct=[0 comp(1,5) ccd(1); ct];
end

if(ccd(end)>ct(end,3))
    ct=[ct; 0 comp(end,6) ccd(end)];
end

offset=ct(:,3)-ct(:,2);
fit=interp1(ct(:,3),offset,ccd);
a=str2num(get(handles.targ_start,'String'));
ind=min(find(ccd-fit>=a)); 
if(~isempty(ind) & ind>1)
    fita=interp1(ccd(ind-1:ind)-fit(ind-1:ind),fit(ind-1:ind),a);
    ccda=a-fita;
else
    fita=fit(1);
    ccda=ccd(1);
end

x1=[ccd(1):1:ccd(end)]';

try
    fit=interp1(ct(:,3),offset,x1);
catch
    disp('Error: Could not interpolate.')
    return
end

fit=[fit; fit(end)];
x1=[x1; 1.2*max(ccd)];



newct=oldct;
if(get(handles.top_fixed,'Value')==1)
    topoff=ccda-top-fita;
    fit=fit+topoff;
    if(~isempty(newct))
        newct(:,2)=newct(:,2)-topoff;
    end
    handles.topoff=topoff;
else
    handles.topoff=0;   
end
ct=[newct;newtie];

ccdoffset=[x1 fit];
H.ccdoffset=ccdoffset;
H.mcdtie=ct;
H.eqn=[];
handles.H=H;
guidata(h,handles);


%---------------------------------------------------------------------------------------
function varargout = poly_fit(h, handles, varargin)
% Find smooth polynomial transformation from ccd to coretop mbsf

H=handles.H;
ccd=H.ccd;
ct=handles.mcdtie_no_offset;

use_eqn=0;
if(~isempty(varargin))
    use_eqn=varargin{1};
end


oldct=ct;
newtie=handles.newtie;
ct=[ct;newtie];

topoff=0;
top=ccd(1);
if(get(handles.top_fixed,'Value')==1)
    top=str2num(get(handles.mcd_start,'String'));
    if(~isempty(top))
        topoff=ccd(1)-top;
    end
end

if(use_eqn==1)
    P=H.eqn';
    mu=H.mu;
    order=max(size(P))-1;
    X=[];
    x1=[min(ccd):2:max(ccd)]';           % x values are in cumulative depth
    x=(x1-mu(1))/mu(2);                % and are centered and scaled to improv fit 
    for i=0:order
        X=[x.^i X];
    end
    fit=P*X';            % y values (estimated offset) for best-fit polynomial
    
    df=fit(2:end)-fit(1:end-1);         % first derivative of estimated offset as function of cumulative depth
    ind1=find(fit'<0 & x1(1:end)<20);   % look for negative offsets in top 20 meters
    sizeneg=size(ind1);
    if(~isempty(ind1))
        % replace negative offstes w/ linear offset from 0 at composite top to first non-neg. offset
        fit(1:ind1(end))=interp1([0 x1(ind1(end)+1)],[0 fit(ind1(end)+1)],x1(1:ind1(end)));
    end
    
    ind2=find(df'<0 & x1(1:end-1)>.85*max(ccd)); % look for decreasing offsets in last 15% of section 
    if(~isempty(ind2))
        fit(ind2(1):end)=fit(ind2(1));     % replace any decreasing offsets with the last non-decreasing offset
    end
    
    
else
    if(isempty(ct))
        H.ccdoffset=[ccd(1) topoff; max(ccd) topoff];
        H.eqn=topoff;
        handles.H=H;
        guidata(h,handles);
        return
    end
    
    mx=max(ct(find(ct(:,1)~=1),3));   % last mcdtie, excluding target core tops
    if(isempty(mx))
        mx=max(ct(:,3));
    end
    x1=[min(ccd):2:mx]';           % x values are in cumulative depth
    ct=ct(find(ct(:,3)<=1.2*mx),:);
    [b,i,j]=unique(ct(:,3));          % indices which exclude any possible duplicates
    sz=max(size(i));
    
    order=get(handles.edit_polyorder,'String');
    offset=ct(i,3)-ct(i,2);
    
    if(~isempty(str2num(get(handles.edit_polyorder,'String'))))
        % find best-fit polynomial of user-defined order 
        order=str2num(get(handles.edit_polyorder,'String'));
        if(order>sz-1)
            order=sz-1;
            set(handles.edit_polyorder,'String',num2str(order))
            disp('Warning: Maximum polynomial order is number of tie points - 1')
        end
        [P,S,mu]=polyfit(ct(i,3),offset-offset(1),order); 
        X=[];
        x=(x1-mu(1))/mu(2);                % and are centered and scaled to improv fit 
        for i=0:order
            X=[x.^i X];
        end
        
    else  % find best-fit polynomials up to order 10
        
        nr=[];
        if(sz<11)
            maxn=sz-1;
        else
            maxn=10;
        end
        for n=1:maxn
            [P,S,mu]=polyfit(ct(i,3),offset-offset(1),n);   
            nr(n)=S.normr;                 % save norm of the residuals for each fit
        end
        if(maxn>1)
            % Find smallest adequate number of terms
            % by selecting the last n which decreases the residuals by at least 1.2% relative to n-1
            n=max(find((nr(1:end-1)-nr(2:end))./nr(1:end-1)>.012))+1;
            if(isempty(n))
                n=maxn;
            end
        else
            n=1;
        end
        x=(x1-mu(1))/mu(2);                % and are centered and scaled to improv fit 
        X10=[x.^10 x.^9 x.^8 x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x x.^0]; 
        [P,S,mu]=polyfit(ct(i,3),offset-offset(1),n);
        X=X10(:,end-n:end);
    end

    fit=P*X';            % y values (estimated offset) for best-fit polynomial
    
    %Fix potential problems at top and bottom of composite, where polynomial less well constrained
    df=fit(2:end)-fit(1:end-1);         % first derivative of estimated offset as function of cumulative depth
    ind1=find(fit'<0 & x1(1:end)<20);   % look for negative offsets in top 20 meters
    sizeneg=size(ind1);
    if(~isempty(ind1))
        % replace negative offstes w/ linear offset from 0 at composite top to first non-neg. offset
        fit(1:ind1(end))=interp1([0 x1(ind1(end)+1)],[0 fit(ind1(end)+1)],x1(1:ind1(end)));
    end
    
    ind2=find(df'<0 & x1(1:end-1)>.85*mx); % look for decreasing offsets in last 15% of section 
    if(~isempty(ind2))
        fit(ind2(1):end)=fit(ind2(1));     % replace any decreasing offsets with the last non-decreasing offset
    end
    
    while(max(diff(fit)>=1)& n>0)    % if offset increases too quickly (.5 m/m), use lower order polynom.
        n=n-1;
        [P,S,mu]=polyfit(ct(i,3),offset-offset(1),n);
        X=X10(:,end-n:end);
        fit=P*X';
        df=fit(2:end)-fit(1:end-1);
        ind1=find(fit'<0 & x1(1:end)<20);
        if(~isempty(ind1))
            fit(1:ind1(end))=interp1([0 x1(ind1(end)+1)],[0 fit(ind1(end)+1)],x1(1:ind1(end)));
        end
        ind2=find(df'<0 & x1(1:end-1)>.85*mx);
        if(~isempty(ind2))
            fit(ind2(1):end)=fit(ind2(1));
        end
        
        df=fit(2:end)-fit(1:end-1);         % first derivative of estimated offset as function of cumulative depth
        ind1=find(fit'<0 & x1(1:end)<20);   % look for negative offsets in top 20 meters
        sizeneg=size(ind1);
        if(~isempty(ind1))
            % replace negative offstes w/ linear offset from 0 at composite top to first non-neg. offset
            fit(1:ind1(end))=interp1([0 x1(ind1(end)+1)],[0 fit(ind1(end)+1)],x1(1:ind1(end)));
        end
        
        ind2=find(df'<0 & x1(1:end-1)>.85*mx); % look for decreasing offsets in last 15% of section 
        if(~isempty(ind2))
            fit(ind2(1):end)=fit(ind2(1));     % replace any decreasing offsets with the last non-decreasing offset
        end
        
    end
    fit=fit+offset(1);
    P(end)=P(end)+offset(1);
end

disp('Equation coefficients:')
P

fit=[fit fit(end)]';       % create constant offset 20% beyond end of ccd
x1=[x1; 1.2*max(ccd)];

a=str2num(get(handles.targ_start,'String'));
fitccd=interp1(x1,fit,ccd);
ind=min(find(ccd-fitccd>=a));
if(~isempty(ind) & ind>1)
    fita=interp1(ccd(ind-1:ind)-fitccd(ind-1:ind),fitccd(ind-1:ind),a);
    ccda=a-fita;
else
    fita=fitccd(1);
    ccda=ccd(1);
end


newct=oldct;
if(get(handles.top_fixed,'Value')==1)
    topoff=ccda-top-fita;
    fit=fit+topoff;
    newct(:,2)=newct(:,2)-topoff;
    handles.topoff=topoff;
else
    handles.topoff=0;   
end

ct=[newct;newtie];
H.mcdtie=ct;

s=['Polynomial fit has order ' num2str(max(size(P))-1)];
disp(s)

H.ccdoffset=[x1 fit];
H.eqn=P;
H.mu=mu;
handles.H=H;
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = push_update_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
if(get(handles.radio_poly,'Value')==1)
    poly_fit(h,handles)
else
    line_interp(h,handles)   
end
plot_offset(h,handles,7);
set(handles.push_update,'Enable','off');


% --------------------------------------------------------------------
function varargout = push_clear_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
H=handles.H;
mx=get(handles.targ_end,'String');
mx=str2num(mx);
if(~isempty(mx))
    H.ccdoffset=[0 0; mx 0];
elseif (~isempty(H.ccd))
    H.ccdoffset=[0 0; max(H.ccd) 0];
else
    H.ccdoffset=[0 0; 500 0];
end
H.mcdtie=[];
handles.mcdtie_no_offset=[];
handles.newtie=[];
handles.H=H;
guidata(h,handles);
plot_offset(h,handles,7);
set(handles.push_update,'Enable','on');


% --------------------------------------------------------------------
function varargout = push_insert_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
%disp('Insert')
%mcdtie= hole, mbsf/mcd, ccd
newpt=insert_mcdtie;
typept=newpt(4);
if(isnan(typept))
    return
end
newpt=newpt(1:3);
H=handles.H;
mcdtie=H.mcdtie;
comp=H.comp;
ccd=H.ccd;

if(typept==0)   % if hole was entered by letter, convert to number
    hole=newpt(1);
    for i=1:H.numholes
        j=eval(['get(H.name' num2str(i) ',''String'')']);
        j=double(j);
        if(hole==j | hole==j+32 | hole==j-32)
            newpt(1)=i;
            typept=1;
            break
        end
    end
    if(typept==0)
        s=['Error: Could not find hole ' char(hole)];
        disp(s)
        return
    end
end

if(typept==1)
    if(newpt(1)==1)
        xm=load(get(H.sig1,'String'));   % original unmatched hole data (mbsf)
        xm(:,3:4)=xm(:,1:2);
    else
        try
            m=['get(H.match' num2str(newpt(1)) ',''String'')'];
            xm=load(eval(m),'-ASCII');
        catch
            s=['Could not find hole ' num2str(newpt(1))];
            disp(s)
            return
        end
    end
 
    mbsf=newpt(2);
    if(mbsf>max(xm(:,2)))
        disp('Error: Depth not found in requested hole.')
        return
    end
    [mx,mc1,c1]=interpc(mbsf,1,[],xm);
    ind=max(find(comp(:,5)<mx));
    %comp(ind,:)
    %ccd(ind:ind+1)
    try
        newpt(2)=newpt(3);
        % Find CCD which corresponds to user requested mbsf
        newpt(3)=interp1(comp(ind,5:6),[ccd(ind) ccd(ind+1)],mx);
    catch
        disp('Error: Could not find equivalent CCD')
        [mbsf mx]
        return
    end
        
else
    newpt=[0 newpt(2) newpt(1)];
end
mcdtie=[mcdtie; newpt];
H.mcdtie=mcdtie;
handles.H=H;

newtie=handles.newtie;
newtie=[newtie; newpt];
handles.newtie=newtie;

guidata(h,handles);
plot_offset(h,handles,7);
set(handles.push_update,'Enable','on');


% --------------------------------------------------------------------
function varargout = top_fixed_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
if(get(handles.top_fixed,'Value')==1)
    set(handles.top_auto,'Value',0);
    f1=get(handles.mcd_start,'String');
    if(isempty(f1))
        set(handles.mcd_start,'String','0');
    end
else
    set(handles.top_auto,'Value',1);
end
set(handles.push_update,'Enable','on');

% --------------------------------------------------------------------
function varargout = top_auto_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
if(get(handles.top_auto,'Value')==1)
    set(handles.top_fixed,'Value',0);
else
    set(handles.top_fixed,'Value',1);
    f1=get(handles.mcd_start,'String');
    if(isempty(f1))
        set(handles.mcd_start,'String','0');
    end
end
set(handles.push_update,'Enable','on');

% --------------------------------------------------------------------
function varargout = mcd_start_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
f1=get(handles.mcd_start,'String');
if(isempty(f1))
    set(handles.top_auto,'Value',1);
    set(handles.top_fixed,'Value',0);
else
    if(~isempty(str2num(f1)))
        set(handles.top_auto,'Value',0);
        set(handles.top_fixed,'Value',1);
    else
        set(handles.top_auto,'Value',1);
        set(handles.top_fixed,'Value',0);
    end
end
set(handles.push_update,'Enable','on');


% --------------------------------------------------------------------
function varargout = targ_start_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
f1=get(handles.targ_start,'String');
if(isempty(f1) | isempty(str2num(f1)))
    comp=handles.H.comp;
    set(handles.targ_start,'String',num2str(comp(1,5)));
end
set(handles.push_update,'Enable','on');
plot_offset(h,handles)

% --------------------------------------------------------------------
function varargout = targ_end_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
f1=get(handles.targ_end,'String');
if(isempty(f1) | isempty(str2num(f1)))
    comp=handles.H.comp;
    set(handles.targ_end,'String',num2str(comp(end,6)));
end
set(handles.push_update,'Enable','on');
plot_offset(h,handles)

% --------------------------------------------------------------------
function varargout = radio_linear_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.
if(get(handles.radio_linear,'Value')==1)
    set(handles.radio_poly,'Value',0);
else
    set(handles.radio_poly,'Value',1);
    if(isempty(get(handles.edit_polyorder,'String')))
        set(handles.edit_polyorder,'String','auto');
    end
end
set(handles.push_update,'Enable','on');

% --------------------------------------------------------------------
function varargout = radio_poly_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.
if(get(handles.radio_poly,'Value')==1)
    set(handles.radio_linear,'Value',0);
    if(isempty(get(handles.edit_polyorder,'String')))
        set(handles.edit_polyorder,'String','auto');
    end
else
    set(handles.radio_linear,'Value',1);
end
set(handles.push_update,'Enable','on');



% --------------------------------------------------------------------
function varargout = edit_polyorder_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
set(handles.push_update,'Enable','on');

% --------------------------------------------------------------------
function varargout = radio_saveeqn_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.

% --------------------------------------------------------------------
function varargout = radio_savefull_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.

% --------------------------------------------------------------------
function varargout = radio_savetie_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.

% --------------------------------------------------------------------
function varargout = save_filename_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.


% --------------------------------------------------------------------
function varargout = push_save_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

fname=get(handles.save_filename,'String');
H=handles.H;
if(isempty(fname))
    disp('Error: Could not save- no filename given')
    set(handles.save_filename,'String','ERROR');
else
    eqn=H.eqn';
    mu=H.mu;  % mean and std. dev. for x from polyfit
    if(get(handles.radio_saveeqn,'Value')==1)
        if isempty(eqn)
            disp('Error: Equation not saved because offset is not a polynomial')
        else
            eqn=flipud(eqn);
            n=[0:max(size(eqn))-1]';
            eqn=[mu(1) mu(2); n eqn];
            save([fname '.eqn'],'eqn','-ASCII');
        end
    end
    if(get(handles.radio_savefull,'Value')==1)
        offset=H.ccdoffset;
        save([fname '.offset'],'offset','-ASCII');
    end
    if(get(handles.radio_savetie,'Value')==1)
        tiepts=H.mcdtie;
        save([fname '.mcdtie'],'tiepts','-ASCII');
    end

end

% --------------------------------------------------------------------
function varargout = push_load_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

%disp('Load')
fname=getconf(8);
if(strcmp(fname,'cancel'))
    return
end

H=handles.H;
[path,name,ext,ver] = fileparts(fname);
if(strcmp(ext,'.eqn'))
    eqn=load(fname);
    H.mu=eqn(1,:);
    eqn=eqn(2:end,:);
    [temp,ind]=sort(eqn(:,1));
    eqn=flipud(eqn(ind,2));
    H.eqn=eqn;
    handles.H=H;
    guidata(h,handles);
    poly_fit(h, handles,1)
    set(handles.radio_linear,'Value',0);
    set(handles.radio_poly,'Value',1);
    set(handles.edit_polyorder,'String',num2str(max(size(eqn))-1));
    set(handles.push_update,'Enable','off');
    
elseif(strcmp(ext,'.mcdtie'))
    mcdtie=load(fname);
    H.mcdtie=mcdtie;
    ind=find(mcdtie(:,1)==0);
    if isempty(ind)
        handles.newtie=[];
    else
        handles.newtie=mcdtie(ind,:);
    end        
    ind=find(mcdtie(:,1)~=0);
    if isempty(ind)
        handles.mcdtie_no_offset=[];
    else
        handles.mcdtie_no_offset=mcdtie(ind,:);
    end        
    set(handles.push_update,'Enable','on');
    handles.H=H;
    guidata(h,handles);

elseif(strcmp(ext,'.offset'))
    H.ccdoffset=load(fname);
    %size(H.ccdoffset)
    set(handles.push_update,'Enable','off');
    handles.H=H;
    guidata(h,handles);
    
else    
    disp('Error: Filename must have the extension .eqn, .mcdtie, or .offset')
    return
end

plot_offset(h, handles,4)

% --------------------------------------------------------------------
function varargout = radio_ccdaxis_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_ccdaxis.
if(get(handles.radio_ccdaxis,'Value')==1)
    set(handles.radio_mcdaxis,'Value',0);
else
    set(handles.radio_mcdaxis,'Value',1);
end
plot_offset(h, handles,4)

% --------------------------------------------------------------------
function varargout = radio_mcdaxis_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_mcdaxis.
if(get(handles.radio_mcdaxis,'Value')==1)
    set(handles.radio_ccdaxis,'Value',0);
else
    set(handles.radio_ccdaxis,'Value',1);
end
plot_offset(h, handles,5)

% --------------------------------------------------------------------
function varargout = push_close_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.push_close.

answer=confirm(6);
if(answer==1)
    handles.apply=0;
    H=handles.H;
    H.test=2;
    handles.H=H;
    guidata(h,handles);
    uiresume(handles.mcdfigure);
end

% --------------------------------------------------------------------
% Return mapping function to Autocomp 
function varargout = push_apply_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

answer=confirm(5);
if(answer==1)
    handles.apply=1;
    H=handles.H;
    H.test=3;
    handles.H=H;
    guidata(h,handles);
    uiresume(handles.mcdfigure);
end

% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
% Stub for ResizeFcn of the figure handles.figure1.


function varargout = autocomp(varargin)
% AUTOCOMP Application M-file for autocomp.fig
%    FIG = AUTOCOMP launch autocomp GUI.
%    AUTOCOMP('callback_name', ...) invoke the named callback.
% tar czfvp Match-2.2.1.tgz Match-2.2.1

% Last Modified by GUIDE v2.0 26-Apr-2007 15:26:05

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    disp(' ')
    disp('Welcome to the Automated Compositing Tool (v2.3)')
    disp('Created by Lorraine Lisiecki, 2003, with funding from')
    disp('a Schlanger Ocean Drilling Fellowship (NSF-USSSP).')
    disp(' ')

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.numholes=1;
    handles.offset=0;
    handles.color=[];
    handles.comp=[];
    handles.cds=[];
    handles.tie=[];
    handles.table=[];
    handles.section=[];
    % the following are set in update or push_editmcd
    handles.apply=[];
    handles.mcdtie=[];
    handles.ccd=[];      % CCD is cumulative composite depth, which is essentially "MCD"
                         % A correction is applied to CCD/MCD to make the CMCD (or NMCD)
                         % depth scale which better fits coretop MBSF measurements 
    handles.mu=[];
    handles.eqn=[];
    handles.ccdoffset=[];
    handles.top=[];
    handles.start=[];
    handles.textx=.95;
    
    %handles.convert_window=0;
    
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
function varargout = sig1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig1.
% --------------------------------------------------------------------
function varargout = sig2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig2.
set(handles.msig2,'String',[get(handles.sig2,'String') '.new']);
% --------------------------------------------------------------------
function varargout = sig3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig3.
set(handles.msig3,'String',[get(handles.sig3,'String') '.new']);
% --------------------------------------------------------------------
function varargout = sig4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.sig4.
set(handles.msig4,'String',[get(handles.sig4,'String') '.new']);

% --------------------------------------------------------------------
function varargout = msig2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.msig2.
% --------------------------------------------------------------------
function varargout = msig3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.msig3.
% --------------------------------------------------------------------
function varargout = msig4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.msig4.

% --------------------------------------------------------------------
function varargout = match2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.match2.
[p,n,e,v]=fileparts(get(handles.match2,'String'));
if(~isempty(n) & isempty(e))
    set(handles.match2,'String',[get(handles.match2,'String') '.match']);
end
% --------------------------------------------------------------------
function varargout = match3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.match3.
[p,n,e,v]=fileparts(get(handles.match3,'String'));
if(~isempty(n) & isempty(e))
    set(handles.match3,'String',[get(handles.match3,'String') '.match']);
end
% --------------------------------------------------------------------
function varargout = match4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.match4.
[p,n,e,v]=fileparts(get(handles.match4,'String'));
if(~isempty(n) & isempty(e))
    set(handles.match4,'String',[get(handles.match4,'String') '.match']);
end

% --------------------------------------------------------------------
function varargout = cdsname_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cdsname.
% --------------------------------------------------------------------
function varargout = name1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.name1.
% --------------------------------------------------------------------
function varargout = name2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.name2.
% --------------------------------------------------------------------
function varargout = name3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.name3.
% --------------------------------------------------------------------
function varargout = name4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.name4.

% --------------------------------------------------------------------
function varargout = plot_height_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
y=ylim(handles.axes1);
yr=str2num(get(handles.plot_height,'String'));
if(~isempty(yr) & yr>0)
    ylim(handles.axes1, [y(1) y(1)+yr]);
    ylim(handles.axes2, [y(1) y(1)+yr]);
end

% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
y=ylim(handles.axes1);
v=get(handles.slider1,'Value');
maxv=get(handles.slider1,'Max');
step=y(1)-(maxv-v);
ylim(handles.axes1, [y(1)-step y(2)-step]);
ylim(handles.axes2, [y(1)-step y(2)-step]);
%xlim(handles.axes1,'auto');



% --------------------------------------------------------------------
% Create a composite section using an automated algorithm, plots results, 
% and calls update_Callback to compute new depth
% --------------------------------------------------------------------
function varargout = create_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.create.

disp('Finding best composite section. This may take ~1 minute.')

if(isempty(get(handles.sig1,'String')))
    disp('Error: No target signal')
    return
else
    %disp('1')
    hole{1}=load(get(handles.sig1,'String'));   % original unmatched hole data (mbsf)
    mhole{1}=hole{1};                           % hole data aligned to target hole
    match{1}=hole{1};                           % conversion from hole to target
    match{1}(:,3:4)=match{1}(:,1:2);
end
if(isempty(get(handles.sig2,'String')))
    disp('Error: No additional signal')
    return
end
s=[];
for i=2:4
    %i
    sig=['get(handles.sig' num2str(i) ',''String'')'];
    if(~isempty(eval(sig)))
        msig=['get(handles.msig' num2str(i) ',''String'')'];
        if(isempty(eval(msig)))
            s=['Error: No matching signal for', get(['handles.sig' num2str(i)],'String')];
            disp(s)
        end
        m=['get(handles.match' num2str(i) ',''String'')'];
        if(isempty(eval(m)))
            s=['Error: No match file for', get(['handles.sig' num2str(i)],'String')];
            disp(s)
        end
        if(isempty(s))
            hole{end+1}=load(eval(sig),'-ASCII');     % original unmatched hole data (mbsf)
            mhole{end+1}=load(eval(msig),'-ASCII');   % hole data aligned to target hole
            match{end+1}=load(eval(m),'-ASCII');      % conversion from hole to target
        end
    end
end
num=max(size(hole));                              % number of holes
step=str2num(get(handles.step,'String'));         % minimum segment size in composite, read from GUI 

d=min(mhole{1}(:,2));     % depth in target hole
c=mhole{1}(1,1);          % target core number
comp=[];
dmax=max(mhole{1}(:,2));  % maximum target depth

while(d+step<=max(mhole{1}(:,2)))    % step down entire length of target hole in segments
    i=find(mhole{1}(:,1)==c);        % find indices of current target core

    % if current target core and start of next core are above top of current 
    % segment (d), find deepest target core which starts at or before segment top.
    if(max(mhole{1}(i,2))<d & mhole{1}(i(end)+1,2)<=d);
        i=max(find(mhole{1}(:,2)<=d));
        c=mhole{1}(i,1);    % deepest core at top of segment
        %new_hole=1;         % badly named, indicates new target core (for tracking target gaps)
    else
        i=max(find(mhole{1}(i,2)<=d));  %index in current target core at or before current segment
        %new_hole=0;  
    end
    
    d2=d+step;      % current segment is d to d2
    
    cind=find(mhole{1}(:,1)==c);
    if(min(mhole{1}(cind,2))>d | max(mhole{1}(cind,2))<d2)
        targ_gap=1;  % target hole contains gap in segment
        %[c d d2 min(mhole{1}(cind,2)) max(mhole{1}(cind,2))]
    else
        targ_gap=0;
    end
    
    i2=min(find(mhole{1}(:,2)>=d2 & mhole{1}(:,1)==c));  % find end of segment in current core
    if(isempty(i2))
        i2=min(find(mhole{1}(:,2)>=d2));    % first index at end of segment
    end
    c2=mhole{1}(i2,1);   %core at end of segment
    new_d=[];
    
    % descriptors of each segment considered for the composite 
    valid=[];     % segment contains no core gapd
    edge=[];      % segment avoids core ends, where distortion is greatest (2m from top, 1m from bottom)
    score=[];     % measure of similarity to other segments
    speed=[];     % rate of extension relative to target
    mhj=[];       % source hole for segment
    mcj=[];       % source core for segment
    mdj=[];       % top of segment, depth in source hole
    md2j=[];      % end of segment, depth in source hole
    my=[];

    % for each hole, find if segment is valid for inclusion in composite and/or near core
    % and assign a score based on proximity to core middle
    %[d d2]
    for j=1:num
        
        h=j;
        di0=find(match{h}(:,4)==d);
        dind=[1:max(size(match{h}(:,4)))]';
        di1=find(match{h}(dind,4)>d & dind>1);
        % find (di1-1,di1) pairs which bracket d and lie in same source core 
        di2=find(match{h}(di1-1,4)<d & match{h}(di1,1)==match{h}(di1-1,1));
        di=[di0; di1(di2)];

        d2i0=find(match{h}(:,4)==d2);
        d2i1=find(match{h}(dind,4)>d2 & dind>1);
        % find (d2i1-1,d2i1) pairs which bracket d2 and lie in same source core 
        d2i2=find(match{h}(d2i1-1,4)<d2 & match{h}(d2i1,1)==match{h}(d2i1-1,1));
        d2i=[d2i0; d2i1(d2i2)];
        
        % for each index di bracketing d, find corresponding index d2i in same core (d & d2 in same core)
        dd2i=0*di;
        if(~isempty(di))
            for i=1:max(size(di))
                ii=find(match{h}(d2i,1)==match{h}(di(i),1));  % all d2i in dame core as di(i)
                if(max(size(ii))>1)   % if multiple soln's find the one that also falls in same target core
                    k=find(match{h}(d2i(ii),3)==match{h}(di(i),3));  
                    ii=ii(k);
                end
                if(max(size(ii))>1)
                    disp('Warning: multiple possible segments for hole, depth:')
                    h
                    match{h}(d2i(ii),:)
                end
                if(~isempty(ii))
                    dd2i(i)=d2i(ii(1));
                end
            end
            if(~isempty(dd2i))
                % set di to only valid solutions, and d2i to matching indices for d2
                di=di(find(dd2i~=0));
                d2i=dd2i(find(dd2i~=0));
            else
                di=[];
                d2i=[];
            end
        end
        
        % set valid, edge, score, speed and location parameters for segment(s)
        if (isempty(di))            % INVALID segments (contain gaps or don't exist)
            valid(end+1)=0;
            edge(end+1)=1;
            score(end+1)=10;
            speed(end+1)=10;
            mhj(end+1)=h;
            % Warning: mcj,mdj, and mdj2 are only approximate for segments containing gaps
            if ~isempty(di)
                mcj(end+1)=match{h}(di(1),1);
            else
                mcj(end+1)=NaN;
                mdj(end+1)=NaN;
            end
            if ~isempty(d2i)
                md2j(end+1)=match{h}(d2i(1),1);
            else
                md2j(end+1)=NaN;
            end
            if(h==1)
                mcj(end)=c;
                mdj(end)=d;
                md2j(end)=d2;
            end
            my(end+1,:)=0*[0:1/20:1]+NaN;               % divide segment into 20 parts
        else                        % VALID segments (don't contain gaps)
            for i=1:max(size(di))
                valid(end+1)=1;
                mhj(end+1)=h;
                mcj(end+1)=match{h}(di(i),1);
                try
                    mdj(end+1)=interp1(match{h}(di(i)-1:di(i),4),match{h}(di(i)-1:di(i),2),d);
                catch
                    mdj(end+1)=match{h}(di(i),2);
                    %disp('Warning: error finding mdj')
                    %[d h]
                    %match{h}(di(i),:)
                end
                try
                    md2j(end+1)=interp1(match{h}(d2i(i)-1:d2i(i),4),match{h}(d2i(i)-1:d2i(i),2),d2);
                catch
                    md2j(end+1)=match{h}(d2i(i),2);
                    %disp('Warning: error finding mdj2')
                    %[d2 h]
                    %match{h}(d2i(i),:)
                end
                md=mdj(end);
                md2=md2j(end);
                
                if(md2<=md)
                    valid(end)=0;
%                     disp('Error: negative segment')
%                     [mhj(end) mcj(end) md md2 d d2]
%                     di
%                     match{h}(di,:)
%                     d2i
%                     match{h}(d2i,:)
                end
                
                % set edge, speed, and score
                speed(end+1)=(md2-md)/(d2-d); % extension in segment relative to target hole
                mcij=find(hole{j}(:,1)==mcj(end));    % indices of core for segment
                x=[md:(md2-md)/20:md2];               % divide segment into 20 parts
                if (isempty(x))
                    x=[md:(md2+1-md)/20:md2+1];               % divide segment into 20 parts
                end
                try
                    % get interpolated signal values in segment for hole j 
                    my(end+1,:)=interp1(hole{j}(mcij,2),hole{j}(mcij,3),x);
                catch
                    my(end+1,:)=0*x+NaN;
                end
                mcmin=min(hole{j}(mcij,2));  % mbsf for core top in hole j
                mcmax=max(hole{j}(mcij,2));  % mbsf for core end in hole j
                if(md-mcmin<2)    % segments most distorted in top 2 meters of core
                    edge(end+1)=-1;   % indicates that segment within 2 meter of core top
                elseif(mcmax-md2 < 1)
                    edge(end+1)=1;    % indicates that segment ends within 1 meter of core end
                else
                    edge(end+1)=0;
                end
                f=(md-mcmin)/(mcmax-mcmin);  % fraction of depth down core
                if edge(end)==1
                    score(end+1)=-1*(f-.5)^2;    % score based on distance from middle of core
                else
                    score(end+1)=(md-mcmin);    % score based on distance from core top
                end
            end
        end
    end
  
    % consider all valid segments that are more than 1 m from core top or end
    choices=find(valid==1 & edge==0);
%     if(max(size(valid))>5)
%         valid 
%         choices
%         speed
%         score
%         my(:,1)
%         size(my)
%     end
    if(isempty(choices))
        % if none found, consider all valid segments
        choices=find(valid==1);
        if(~isempty(choices))
            % select valid choice with minimum distance from center of core
            %choice=choices(find(score(choices)==min(score(choices))));
            choice=choices(find(score(choices)==max(score(choices))));
        else
            % if no valid choices, use target hole
            choice=1;
        end

    % if more than one valid, non-edge choice    
    elseif(max(size(choices))>2)
        choices=find(valid==1);
        
        % if multiple segments have similar rates of extension, use those holes
        % i.e., try to eliminate segments which don't match the others
        spch=find(speed(choices)<1.5 & speed(choices)>.67);
        if(max(size(spch))<2)
            spch=find(speed(choices)>=1.5);
            if(max(size(spch))<2)
                spch=find(speed(choices)<=.67);
                if(max(size(spch))<2)
                    spch=find(speed(choices)<1.5 & speed(choices)>.67);
                    % if no 2 holes alike, revert back to 0.67-1.5
                end  
            end  
        end  
        %spch
        choices=choices(spch);

        % if more than 2 segments have similar extensions, calculate correlation
        % coefficients of segment signal for each pair of holes and choose the
        % 2 with the best correlation
        if(max(size(choices))>2)
            cmax=-1;
            for j=1:max(size(choices))-1
                for k=j+1:max(size(choices))
                    temp=corrcoef(my(choices(j),:),my(choices(k),:));
                    if(temp(1,2)>cmax)
                        cmax=temp(1,2);
                        cc=[j k];
                    end
                end
            end
            choices=choices(cc);
        end
        
        % among remaining valid choices, choose the one farthest from core top  
        if(~isempty(choices))
            %choice=find(score==min(score(choices)));
            choice=find(score==max(score(choices)));
        else  % this error would indicate a bug in selection code above
            disp('Error in selecting choice')
        end
    else   % if 1 or 2 valid, non-edge choices, choose the one farthest from core top 
        choice=find(score==max(score(choices)));
    end
    
    % this should never happen because choice set to target above if no valid choices
    if(isempty(choice))
        disp('Error: No valid choice')
    end

    if(choice(1)<0 | choice(1)>max(size(mhj)))
        disp('Error: cannot find choice')
        choice
    end
    
    % store hole values for selected segment
    comp(end+1,1)=mhj(choice(1));     % hole 
    comp(end,2)=mcj(choice(1));   % core
    comp(end,3)=mdj(choice(1));  % top of segment in hole mbsf
    comp(end,4)=md2j(choice(1)); % end of segment in hole mbsf
    comp(end,5)=d;               % top of segment in target mbsf
    comp(end,6)=d2;              % end of segment in target mbsf
    comp(end,7)=targ_gap;        % 1 if previous segment contains a target gap. UNNECESSARY?? ELIMINATE????

    
    d=d+step;
end

if(comp(1,1)~=1)  % if first segment not from target hole
    comp(1,3)=min(hole{comp(1,1)}(:,2));  %top depth set to hole's first core top
    comp(1,5)=comp(1,3);            %assume hole depth equivalent to same depth in target
end

%comp
%disp('Finished comp')

for j=1:max(size(comp(:,1)))
    if(comp(j,3)>comp(j,4))
        disp('Error: Negative segment in composite')
        comp(j-1:j+1,:)
    end
end

%disp('Plot')

%-----------------------------------------------------
% PLOT COMPOSITE RESULTS
axes(handles.axes1)
cla
% customize button down function inside this plot
set(gca,'DefaulttextButtonDownFcn','autocomp(''text_ButtondownFcn'',gcbo,[],guidata(gcbo))');

offset=3*std(hole{1}(:,3));
meanx=mean(hole{num}(:,3));
%offset=2.7*std(hole{1}(:,3));

for j=1:num
    if(round(j/2)==j/2)
        gr=.5;
    else
        gr=0;
    end
    color(j,:)=[0 gr (j-1)/(num-1)];
    c=mhole{j};
    plot(c(:,3)+(j-1)*offset,c(:,2),'Color',color(j,:));
    hold on;
    pc=match{j}(1,1);
    for k=1:max(size(match{j}(:,1)))
        if(match{j}(k,1)~=pc)
            ind=max(find(mhole{j}(:,2)<=match{j}(k,4)));
            plot(mhole{j}(ind,3)+(j-1)*offset, mhole{j}(ind,2),'o','Color',color(j,:));
            pc=match{j}(k,1);
        end
    end
end

% customize button down function inside this plot, doing it again just in case....
set(gca,'DefaulttextButtonDownFcn','autocomp(''text_ButtondownFcn'',gcbo,[],guidata(gcbo))');
j=1;
cnt=1;
section=[];
while(j<=max(size(comp)))
    for k=j:max(size(comp))
        if(comp(j,1)~=comp(k,1))
            break;
        end 
    end
    if(k~=max(size(comp)))
        k=k-1;
    end
    h=comp(j,1);
    ind=find(mhole{h}(:,2)>=comp(j,5) & mhole{h}(:,2)<=comp(k,6));
    section(end+1,:)=[cnt j k];     %segment number, first and last lines in comp of segments from same hole
    plot(mhole{h}(ind,3)+(h-1)*offset,mhole{h}(ind,2),'r');
    plot([mean(hole{1}(:,3))-1.25*offset (num+.25)*offset+meanx],...
        [mhole{h}(ind(1),2) mhole{h}(ind(1),2)],'k:');
    text(num*offset+meanx,mhole{h}(ind(1),2),num2str(cnt),'FontSize',7,'VerticalAlignment','top');
    %text((num+1.85)*offset+meanx,mhole{h}(ind(1),2),num2str(cnt),'FontSize',7,'VerticalAlignment','top');
    cnt=cnt+1;
    j=k+1;
end
handles.numholes=num;
handles.offset=offset;
handles.color=color;
handles.comp=comp;
handles.section=section;
handles.textx=num*offset+meanx;
guidata(gcbo, handles);


ymax=str2num(get(handles.plot_height,'String'));
if(isempty(ymax) | ymax<=0)
    ymax=40;
    set(handles.plot_height,'String','40')
end
%axis tight
ylim([0 ymax])
xlim([mean(hole{1}(:,3))-1.25*offset (num+.25)*offset+meanx])
maxv=max(hole{1}(:,2));
set(handles.slider1,'min',0,'max',maxv,'Value',maxv,'sliderstep',[0.04 0.1]);

axis ij;

%disp('UPDATE')
update_Callback(h, eventdata, handles, 0)
disp('Composite section finished.')


% --------------------------------------------------------------------
function varargout = push_editmcd_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radio_editmcd.

handles=guidata(handles.figure1);
comphandles=ccd_mcd(handles);
%test=comphandles.test
if(~isempty(comphandles.apply) & comphandles.apply==1)
    handles=comphandles;
    guidata(gcbo, handles);
    update_Callback(h, eventdata, handles, 1)
end

% --------------------------------------------------------------------------------------
% Stretches/compresses composite to fit coretop mbsf values
% and updates composite section to include user changes
% ---------------------------------------------------------------------------------------
function varargout = update_Callback(h, eventdata, handles, varargin)

load_offset=0;
if(~isempty(varargin))
    if(isempty(varargin{1}))
        load_offset=0;
    else
        load_offset=varargin{1};
    end
end

handles=guidata(handles.figure1);
%comp: 1-hole, 2-core, 3-mbsf start, 4-mbsf stop, 5-target start, 6-target stop, 
%      7-flag for contains target coretop
comp=handles.comp;
size_comp=size(comp);
name=get(handles.cdsname,'String');
%save([name '.notie.comp'],'comp','-ASCII');

hole{1}=load(get(handles.sig1,'String'));
mhole{1}=hole{1};
match{1}=hole{1};
match{1}(:,3:4)=match{1}(:,1:2);
for i=2:4
    sig=['get(handles.sig' num2str(i) ',''String'')'];
    if(isempty(eval(sig)))
        break;
    end
    hole{end+1}=load(eval(sig),'-ASCII');
    msig=['get(handles.msig' num2str(i) ',''String'')'];
    mhole{end+1}=load(eval(msig),'-ASCII');
    m=['get(handles.match' num2str(i) ',''String'')'];
    match{end+1}=load(eval(m),'-ASCII');
end
num=max(size(hole));

if(load_offset==1)
    %ctop=handles.mcdtie;
    ccdoffset=handles.ccdoffset;
    x1=ccdoffset(:,1);
    fit=ccdoffset(:,2)';
    ccd=handles.ccd;
else
    ccd=[];
    handles.mcdtie=[];
end

if(isempty(ccd))
    % calculate cumulative composite depth
    ccd=comp(1,3);
    for i=2:max(size(comp))+1
        ccd(i)=ccd(i-1)+comp(i-1,4)-comp(i-1,3);
    end
end
sizeccd=size(ccd);

if(isempty(handles.mcdtie))
    xm=hole{1};
    ctop=[];   % record of all coretop depths, in original mbsf and cumulative depth 
    % loop through all target cores to create ctop for target
    for c=1:max(xm(:,1))
        ind=find(xm(:,1)==c);       % indices for current core
        if(~isempty(ind))
            mbsf=xm(ind(1),2);        % coretop mbsf depth
            jnd=find(comp(:,6)>mbsf); % comp entries containing and following core top
            if(~isempty(jnd))
                k=jnd(1);               % entry spanning target core top
                df=comp(k,4)-comp(k,3); % length of segment in target mbsf
                %convert coretop depth to  cumulative composite depth
                yc=interp1([comp(k,5) comp(k,6)],[ccd(k) ccd(k)+df],mbsf);
                ctop(end+1,:)=[1 mbsf yc];
            end     
        end
    end
    % find coretops for non-target holes
    for j=2:num
        xm=match{j};
        for c=1:max(xm(:,1))
            ind=find(xm(:,1)==c);
            if(~isempty(ind))
                mbsf=xm(ind(1),2);        % coretop mbsf depth (this hole)
                y=xm(ind(1),4);           % coretop mbsf depth (target)
                jnd=find(comp(:,6)>y);    % comp entries containing and following core top
                if(~isempty(jnd))
                    k=jnd(1);
                    df=comp(k,4)-comp(k,3); % length of segment in hole mbsf
                    %convert coretop depth to  cumulative composite depth
                    yc=interp1([comp(k,5) comp(k,6)],[ccd(k) ccd(k)+df],y);
                    ctop(end+1,:)=[j mbsf yc];
                end     
            end
        end
    end
    handles.mcdtie=ctop;    
end

if(load_offset~=1)    
    %---------------------------------------------------------------------------------------
    % Find smooth polynomial transformation from ccd to coretop mbsf
    mx=max(ctop(find(ctop(:,1)~=1),3));     % last non-target core top
    ct=ctop(find(ctop(:,3)<=1.2*mx),:);     % all coretops where holes overlap
    ct=[1 0 0; 1 .01 .01; 1 .02 .02; ct];   % add extra constraints so fit will be ~0 at hole tops
    [b,i,j]=unique(ct(:,3));                % indices which exclude any possible duplicates
    %ctop(i,:)
    nr=[];
    if(max(size(ct(i,3)))<11)
        maxn=max(size(ct(i,3)))-1;
    else
        maxn=10;
    end
    for n=1:maxn
        [P,S,mu]=polyfit(ct(i,3),ct(i,3)-ct(i,2),n);   % find best-fit polynomials up to order 10
        nr(n)=S.normr;                                 % save norm of the residuals for each fit
    end
    % Find smallest adequate number of terms
    % by selecting the last n which decreases the residuals by at least 1.2% relative to n-1
    n=max(find((nr(1:end-1)-nr(2:end))./nr(1:end-1)>.012))+1;
    x1=[0:2:max(ct(i,3))]';                 % x values are in cumulative depth
    x=(x1-mu(1))/mu(2);                     % and are centered and scaled to improv fit 
    X10=[x.^10 x.^9 x.^8 x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x x.^0]; 
    [P,S,mu]=polyfit(ct(i,3),ct(i,3)-ct(i,2),n);
    X=X10(:,end-n:end);
    fit=P*X';                         % y values (estimated offset) for best-fit polynomial
    %size(fit)
    
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
        [P,S,mu]=polyfit(ct(i,3),ct(i,3)-ct(i,2),n);
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
    end
    
    fit=[fit fit(end)];       % create constant offset 20% beyond end of ccd
    x1=[x1; 1.2*max(ccd)];
    ccdoffset=[x1 fit'];
    handles.mcdtie=ctop;
    handles.eqn=P;
    handles.mu=mu;
end


%rebnd=[x1 fit'];     % composite depth, ccd-mbsf offset
%rebnd=[x1-fit' fit'];     % composite depth, ccd-mbsf offset
%save([name '.adjust'],'rebnd','-ASCII');

%[ccd(2:end)' comp(:,6)]
%[x1(1:3) fit(1:3)']
%[x1(end-2:end) fit(end-2:end)']

%---------------------------------------------------------------------------
% Convert from ccd (d) to adjusted composite depth (dc) in composite section
comp2=comp;
y=[];
d=[];
dc=[];
emptyi=[];
j=1;
while(j<=max(size(comp)))
    % Find all adjacent segments from a single hole 
    for k=j:max(size(comp))
        if(comp(j,1)~=comp(k,1))
            break;
        end 
    end
    if((k~=max(size(comp)) | comp(j,1)~=comp(end,1)) & k>1)
        k=k-1;
    end
    
    h=comp(j,1);
    
    % indices for segment in original hole
    % in case of overlap, take points from upper core
    hind=[];
    cend=hole{h}(1,2)-.1;   
    for ci=comp(j,2):comp(k,2)+1
        ind=find(hole{h}(:,2)>=min(comp(j:k,3)) & ...
            hole{h}(:,2)<=max(comp(j:k,4)) & hole{h}(:,1)==ci);
        if(isempty(hind))
            hind=ind;
            cend=max(hole{h}(hind,2));
        else
            ind=ind(find(hole{h}(ind,2)>cend));
            if(~isempty(ind))
                hind(end+1:end+max(size(ind)))=ind;
                cend=max(hole{h}(hind,2));
            end
        end
    end
       
    if(~isempty(hind))
        md=comp2(j,3);
        md2=comp2(k,4);
        offset1=interp1(x1,fit',ccd(j));
        offset2=interp1(x1,fit',ccd(k+1));
        try
            newd=interp1([md md2],[ccd(j)-offset1 ccd(k+1)-offset2],hole{h}(hind,2));
        catch
            [md, md2, ccd(j), offset1, ccd(k+1), offset2, hole{h}(hind(1),2)]
        end
        dc(end+1:end+max(size(hind)))=newd;                   % composite depth
        y(end+1:end+max(size(hind)))=hole{h}(hind,3);
        comp2(j:k,5)=interp1([comp(j,5) comp(k,6)], ...
            [ccd(j)-offset1 ccd(k+1)-offset2],comp2(j:k,5));
        comp2(j:k,6)=interp1([comp(j,5) comp(k,6)], ...
            [ccd(j)-offset1 ccd(k+1)-offset2],comp2(j:k,6));
        if(isnan(comp2(j,5)))
            disp('ERROR: Could not calculate MCD')
            [comp(j,5) comp(k,6), ccd(j), offset1, ccd(k+1), offset2, comp2(j,5:6)]
        end

    else  % if no data for this segment found in hole, scale mcd based on length of missing material
        disp('Warning: Empty segment in composite section for line numbers-')
        md=comp2(j,3);
        md2=comp2(k,4);
        offset1=interp1(x1,fit',ccd(j));
        offset2=interp1(x1,fit',ccd(k+1));
        comp2(j:k,5)=interp1([d d2],[ccd(j)-offset1 ccd(k+1)-offset2],comp2(j:k,5));
        comp2(j:k,6)=interp1([d d2],[ccd(j)-offset1 ccd(k+1)-offset2],comp2(j:k,6));     
    end
    j=k+1;
end

% if (~isempty(emptyi))
%     emptyi
%     comp(emptyi(1)-1:emptyi(end)+1,:)
%     comp2(emptyi(1)-1:emptyi(end)+1,:)
% end

%save('testcomp','comp','-ASCII');


% Final composite depth section
cds=[dc' y'];
%save([name '.cds'],'cds','-ASCII');

% Concise description of composite section for output upon saving
table=[];
j=1;
while(j<=max(size(comp)))
    for k=j:max(size(comp))
        if(comp(j,1)~=comp(k,1) | comp(j,2)~=comp(k,2))
            break;
        end 
    end
    if(k~=max(size(comp)) & k>1)
        k=k-1;
    end
    table(end+1,:)=[comp(j,1:3) comp(k,4) comp(j,5) comp2(j,5)];
    j=k+1;
end
j=k;
table(end+1,:)=[comp(k,1:2) comp(k,4) comp(k,4) comp(k,6) comp2(j,6)];

% for i=1:max(size(table(:,1)))
%     ind=find(hole{table(i,1)}(:,1)==table(i,2));
%     if(max(hole{table(i,1)}(ind,2))<table(i,4))
%         table(i,4)=max(hole{table(i,1)}(ind,2));
%         %table(i,:)
%     end
% end

%table
%comp(1:10,:)
%comp2(1:10,:)

handles.ccdoffset=ccdoffset;
handles.ccd=ccd;

handles.table=table;
handles.cds=cds;
guidata(gcbo, handles);

axes(handles.axes2);
cla
plot(cds(:,2),cds(:,1));
ylim(ylim(handles.axes1))
axis ij
hold on


% --------------------------------------------------------------------
function varargout = step_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.load.


% --------------------------------------------------------------------
function varargout = load_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.load.
name=get(handles.cdsname,'String');
comp=load([name '.comp']);
a=textread([name '.files'],'%s');
n=find(strcmp(a,'step'));
if(isempty(n))
    disp('Error: could not read file')
    return
end
i=1;
h=0;
label=0;
while(i<n)
    h=h+1;
    if(max(size(a{i}))==1)
        sig=['set(handles.name' num2str(h) ',''String'',a{i})'];
        eval(sig);
        i=i+1;
        label=1;
    end
    s=['handles.sig' num2str(h)];
    set(eval(s),'String',a{i});
    hole{h}=load(a{i},'-ASCII');
    i=i+1;
    if(h>1)
        s=['handles.msig' num2str(h)];
        set(eval(s),'String',a{i});
        mhole{h}=load(a{i},'-ASCII');
        i=i+1;
        s=['handles.match' num2str(h)];
        set(eval(s),'String',a{i});
        match{h}=load(a{i},'-ASCII');
        i=i+1;
    else
        mhole{1}=hole{1};
        match{1}=[hole{1}(:,1:2) hole{1}(:,1:2)];
    end
end
num=h;
if(num<4)
    for i=num+1:4
        s=['handles.name' num2str(i)];
        set(eval(s),'String','');
        s=['handles.sig' num2str(i)];
        set(eval(s),'String','');
        s=['handles.msig' num2str(i)];
        set(eval(s),'String','');
        s=['handles.match' num2str(i)];
        set(eval(s),'String','');
    end
end
set(handles.step,'String',a{n+1});
n=find(strcmp(a,'equation'));
if(~isempty(n))
    fname=a{n+1};
    eqn=load(fname);
    handles.mu=eqn(1,:);
    eqn=eqn(2:end,:);
    [temp,ind]=sort(eqn(:,1));
    eqn=flipud(eqn(ind,2));
    handles.eqn=eqn;
else
    handles.eqn=[];    
end
n=find(strcmp(a,'offset'));
if(~isempty(n))
    fname=a{n+1};
    handles.ccdoffset=load(fname);
    handles.ccd=[];
else
    handles.ccdoffset=[];    
    handles.ccd=[];
end

axes(handles.axes1)
cla
set(gca,'DefaulttextButtonDownFcn','autocomp(''text_ButtondownFcn'',gcbo,[],guidata(gcbo))');

offset=3*std(hole{1}(:,3));
meanx=mean(hole{num}(:,3));

for j=1:num
    if(round(j/2)==j/2)
        gr=.5;
    else
        gr=0;
    end
    color(j,:)=[0 gr (j-1)/(num-1)];
    c=mhole{j};
    plot(c(:,3)+(j-1)*offset,c(:,2),'Color',color(j,:));
    hold on;
    pc=match{j}(1,1);
    for k=1:max(size(match{j}(:,1)))
        if(match{j}(k,1)~=pc)
            ind=max(find(mhole{j}(:,2)<=match{j}(k,4)));
            plot(mhole{j}(ind,3)+(j-1)*offset, mhole{j}(ind,2),'o','Color',color(j,:));
            pc=match{j}(k,1);
        end
    end
end

set(gca,'DefaulttextButtonDownFcn','autocomp(''text_ButtondownFcn'',gcbo,[],guidata(gcbo))');
j=1;
cnt=1;
section=[];
while(j<=max(size(comp)))
  for k=j:max(size(comp))
    if(comp(j,1)~=comp(k,1))
      break;
    end 
  end
  if(k~=max(size(comp)))
    k=k-1;
  end
  h=comp(j,1);
  ind=find(mhole{h}(:,2)>=comp(j,5) & mhole{h}(:,2)<=comp(k,6));
  section(end+1,:)=[cnt j k];
  plot(mhole{h}(ind,3)+(h-1)*offset,mhole{h}(ind,2),'r');
  plot([mean(hole{1}(:,3))-1.25*offset (num+.25)*offset+meanx],...
      [mhole{h}(ind(1),2) mhole{h}(ind(1),2)],'k:');
  text(num*offset+meanx,mhole{h}(ind(1),2),num2str(cnt),'FontSize',7,'VerticalAlignment','top');
  %text((num+1.85)*offset+meanx,mhole{h}(ind(1),2),num2str(cnt),'FontSize',7,'VerticalAlignment','top');
  cnt=cnt+1;
  j=k+1;
end
handles.numholes=num;
handles.offset=offset;
handles.color=color;
handles.comp=comp;
handles.section=section;
handles.textx=num*offset+meanx;
guidata(gcbo, handles);


ymax=str2num(get(handles.plot_height,'String'));
if(isempty(ymax) | ymax<=0)
    ymax=40;
    set(handles.plot_height,'String','40')
end
%axis tight
ylim([0 ymax])
xlim([mean(hole{1}(:,3))-1.25*offset (num+.25)*offset+meanx])
axis ij;
maxv=max(hole{1}(:,2));
set(handles.slider1,'min',0,'max',maxv,'Value',maxv,'sliderstep',[0.04 0.1]);

if(~isempty(handles.ccdoffset))
    %disp('update 1')
    update_Callback(h, eventdata, handles, 1)
else
    update_Callback(h, eventdata, handles, varargin)
end

% --------------------------------------------------------------------
function varargout = save_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.save.
hole{1}=load(get(handles.sig1,'String'));
mhole{1}=hole{1};
match{1}=hole{1};
match{1}(:,3:4)=match{1}(:,1:2);
for i=2:4
    sig=['get(handles.sig' num2str(i) ',''String'')'];
    if(isempty(eval(sig)))
        break;
    end
    hole{end+1}=load(eval(sig),'-ASCII');
    msig=['get(handles.msig' num2str(i) ',''String'')'];
    mhole{end+1}=load(eval(msig),'-ASCII');
    m=['get(handles.match' num2str(i) ',''String'')'];
    match{end+1}=load(eval(m),'-ASCII');
end
num=max(size(hole));

comp=handles.comp;
name=get(handles.cdsname,'String');
save([name '.comp'],'comp','-ASCII');

fc=[name,'.files'];
fid=fopen(fc,'w');

c=get(handles.name1,'String');
if(~isempty(c))
    fprintf(fid,'%c\n',c(1));
    h(1)=c;
else
    h(1)=num2str(i);
end
fprintf(fid,'%s\n',get(handles.sig1,'String'));
signal{1}=get(handles.sig1,'String');
for i=2:4
    sig=['get(handles.name' num2str(i) ',''String'')'];
    c=eval(sig);
    if(~isempty(c))
        fprintf(fid,'%c\n',c(1));
        h(i)=c;
    else
        h(i)=num2str(i);
    end
    sig=['get(handles.sig' num2str(i) ',''String'')'];
    fprintf(fid,'%s\n',eval(sig));
    signal{i}=eval(sig);
    sig=['get(handles.msig' num2str(i) ',''String'')'];
    fprintf(fid,'%s\n',eval(sig));
    sig=['get(handles.match' num2str(i) ',''String'')'];
    fprintf(fid,'%s\n',eval(sig));
end

fprintf(fid,'step %s\n',get(handles.step,'String'));

if(~isempty(handles.eqn))
    fprintf(fid,'equation %s\n',[name '.eqn']);
end
if(~isempty(handles.ccdoffset))
    fprintf(fid,'offset %s\n',[name '.offset']);
end
fclose(fid);

if(~isempty(handles.eqn))
    eqn=handles.eqn;
    size_eqn=size(eqn);
    if(size_eqn(:,1)<size_eqn(:,2))
        eqn=eqn';
    end
    mu=handles.mu;  % mean and std. dev. for x from polyfit
    eqn=flipud(eqn);
    n=[0:max(size(eqn))-1]';
    eqn=[mu(1) mu(2); n eqn];
    %eqn=[n eqn];
    save([name '.eqn'],'eqn','-ASCII');
end
if(~isempty(handles.ccdoffset))
    ccdoffset=handles.ccdoffset;
    save([name '.offset'],'ccdoffset','-ASCII');
end

table=handles.table;

ind=find(table(:,4)<table(:,3));
if(~isempty(ind))
    disp('Suspicious splice removed:')
    table(ind,:)
    table(ind,:)=[];
end

if(~isempty(find(table(:,5)==999 | table(:,6)==999)))
    disp('Warning: Table for composite section construction contains invalid values denoted 999.')
end

fc=[name,'.section'];
fid=fopen(fc,'w');
fprintf(fid,'CMCD = MCD');
for j=2:max(size(eqn))
    fprintf(fid,' + %8.3f x^%2d',-1*eqn(j,2), eqn(j,1));
end
fprintf(fid,'\n\n');
fprintf(fid,'Hole Core Start(MBSF) End(MBSF) Start(MCD) Start(CMCD)\n');
for j=1:max(size(table))
    fprintf(fid,'%4c %4d %11.4f %9.4f %11.4f %9.4f\n',h(table(j,1)), table(j,2:6));
end
fclose(fid);

cds=handles.cds;
save([name '.cds'],'cds','-ASCII');

% This code used to generate tie points for alignments to the composite section, but it is no longer supported
% for i=2:num
%     tp{i}=[];
%     for c=match{i}(1,1):match{i}(end,1)
%         ind=find(match{i}(:,1)==c);
%         ind2=find(hole{i}(:,1)==c);
%         if(isempty(ind))
%             continue
%         end
%         jnd=min(find(newh(:,1)==match{i}(ind(1),3) & newh(:,2)>=match{i}(ind(1),4)));
%         if(~isempty(jnd) & jnd~=1)
%             if(newh(jnd-1,1)==newh(jnd,1))
%                 if(match{i}(ind(1),2)>=min(hole{i}(ind2,2)))
%                     newx=interp1(newh(jnd-1:jnd,2),newh(jnd-1:jnd,4),match{i}(ind(1),4));
%                     if(~isnan(newx))
%                         tp{i}(end+1,:)=[match{i}(ind(1),1:2) 0 newx];
%                     end
%                 end
%             end
%         end
%         for j=1:3
%             k=round(j/3*max(size(ind)));
%             if(k>1 & match{i}(ind(k),2)>max(hole{i}(find(hole{i}(:,1)==c))))
%                 k=k-1;
%             end
%             jnd=min(find(newh(:,1)==match{i}(ind(k),3) & newh(:,2)>=match{i}(ind(k),4)));
%             if(~isempty(jnd) & jnd~=1)
%                 if(newh(jnd-1,1)==newh(jnd,1))
%                     newx=interp1(newh(jnd-1:jnd,2),newh(jnd-1:jnd,4),match{i}(ind(k),4));
%                     if(~isnan(newx))
%                         tp{i}(end+1,:)=[match{i}(ind(k),1:2) 0 newx];
%                     end
%                 end
%             end
%         end
%     end
%     t=tp{i};
%     [p,n,e,v]=fileparts(signal{i});
    %save([n '.mcd.tie'],'t','-ASCII');
    %end


% --------------------------------------------------------------------
function varargout = edit_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit.

% --------------------------------------------------------------------
% Quit program (graphical user interface, autocomp)
function varargout = exit_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit.
delete(handles.figure1);

% --------------------------------------------------------------------
function varargout = delete_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.delete.
if(get(handles.edit,'Value')==1)
    n=str2num(get(handles.delete,'String'));
    if(isempty(n))
        return
    end
    if(n==1)
        disp('Cannot delete first section')
        return
    end
    comp=handles.comp;
    section=handles.section;
    %size(section)
    ni=find(section(:,1)==n);
    if(isempty(ni))
        s=['Error: Section ' num2str(n) ' does not exist.'];
        disp(s);
        return
    end
    a=section(ni,2);
    b=section(ni,3);
    %comp(a-1:b+1,:)
    %section(ni-1:ni+2,:)
    ind=[a:b];
    i=comp(a-1,1);
    oldi=comp(a,1);
    targ=load(get(handles.sig1,'String'));
    if(i==1)
        %hole=targ;
        mhole=targ;
        match=[targ(:,1:2) targ(:,1:2)];
    else
        sig=['get(handles.sig' num2str(i) ',''String'')'];
        msig=['get(handles.msig' num2str(i) ',''String'')'];
        m=['get(handles.match' num2str(i) ',''String'')'];
        %hole=load(eval(sig),'-ASCII');
        mhole=load(eval(msig),'-ASCII');
        match=load(eval(m),'-ASCII');
    end
    if(oldi==1)
        omhole=targ;
    else
        msig=['get(handles.msig' num2str(oldi) ',''String'')'];
        omhole=load(eval(msig),'-ASCII');
    end

    newline=[];
    if(~isempty(ind))
        for j=1:max(size(ind))
            d=comp(ind(j),5);           
            d2=comp(ind(j),6);
            k=min(find(targ(:,2)>=d));
            c1=targ(k,1);
            c2=c1+1;
            %Note: hole should not contain gap.
            mci=find(match(:,3)>=c1-1 & match(:,3)<=c2+2);
            mi1_less=max(find(match(mci,4)<=d & match(mci,3)==c1));
            if(isempty(mi1_less))
                mi1_less=max(find(match(mci,4)<=d));
            end
            mc=match(mci(mi1_less),1);
            mi1_more=min(find(match(mci,4)>d & match(mci,1)==mc));
            mi2_less=max(find(match(mci,4)<=d2 & match(mci,1)==mc));
            mi2_more=min(find(match(mci,4)>d2 & match(mci,1)==mc));
            if(isempty(mi1_more)| mi1_more<=mi1_less)
                mi1_more=min(find(match(mci,4)>d & mci>mi1_less));
            end
            if(isempty(mi2_less))
                mi2_less=max(find(match(mci,4)<=d2));
            end
            if(isempty(mi2_more) | mi2_more<=mi2_less)
                mi2_more=min(find(match(mci,4)>=d2 & mci>mi2_less));               
            end
            %[c1 c2]
            %[mi1_less mi1_more mi2_less mi2_more]
            %mci(mi1_less)
            sub=match(mci(mi1_less:mi1_more),:);
            md=interp1(sub(:,4),sub(:,2),d);
            sub=match(mci(mi2_less:mi2_more),:);
            md2=interp1(sub(:,4),sub(:,2),d2);
            comp(ind(j),1:6)=[i mc md md2 d d2];
        end
    end
    %section(ni-1:ni+1,:)
    section(ni-1,3)=section(ni,3);
    section(ni,:)=[];
    %comp(a-1:b+1,:)
    %disp('deleted:')
    %section(ni-1:end,:)
    handles.comp=comp;
    handles.section=section;
    guidata(gcbo, handles);
   
    hobj1=findobj('String',num2str(n),'Parent',handles.axes1,'Color','k');
    x=get(hobj1,'Position');
    delete(hobj1);
    delete(findobj('YData',[x(2) x(2)],'Parent',handles.axes1));
    axes(handles.axes1)
    hold on;
    ind=find(omhole(:,2)>=comp(a,5) & omhole(:,2)<=comp(b,6));
    color=handles.color;
    plot(omhole(ind,3)+(oldi-1)*handles.offset,omhole(ind,2),'Color',color(oldi,:));
    ind=find(mhole(:,2)>=comp(a,5) & mhole(:,2)<=comp(b,6));
    plot(mhole(ind,3)+(i-1)*handles.offset,mhole(ind,2),'r');

end

% --------------------------------------------------------------------
function varargout = insert_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.insert.
if(get(handles.edit,'Value')==1)
    n=str2num(get(handles.insert,'String'));
    if(isempty(n))
        return
    end
    comp=handles.comp;
    section=handles.section;
    ni=find(section(:,1)==n);
    if(isempty(ni))
        s=['Error: Section ' num2str(n) ' does not exist.'];
        disp(s);
        return
    end
    a=section(ni,2);
    b=section(ni,3);
    step=str2num(get(handles.step,'String'));
    if(comp(b,6)-comp(a,5)<=step)
        s=['Error: Section ' num2str(n) ' is too small to subdivide.'];
        disp(s);
        disp('Decrease step size to allow smaller segments.')
        return
    end

    %comp(a-1:b+1,:)
    %section(ni-1:ni+2,:)
    d=comp(a,5);         
    d2=comp(a,5)+step;
    if(comp(a,6)~=d2)
        i=comp(a,1);
        targ=load(get(handles.sig1,'String'));
        if(i==1)
            mhole=targ;
            match=[targ(:,1:2) targ(:,1:2)];
        else
            msig=['get(handles.msig' num2str(i) ',''String'')'];
            m=['get(handles.match' num2str(i) ',''String'')'];
            mhole=load(eval(msig),'-ASCII');
            match=load(eval(m),'-ASCII');
        end
        k=max(find(targ(:,2)<=d2));
        c2=comp(a,2);
        c1=c2;
        mci=find(match(:,1)==c2+1 | match(:,1)==c2);
        mi2_less=max(find(match(mci,4)<=d2 & match(mci,1)==c2));
        mi2_more=min(find(match(mci,4)>d2 & match(mci,1)==c2));
        new=0;
        if(isempty(mi2_less))
            mi2_less=max(find(match(mci,4)<=d2));
            c2=match(mci(mi2_less),1);
            new=1;
        end
        if(isempty(mi2_more) | mi2_more<=mi2_less)
            mi2_more=min(find(match(mci,4)>=d2 & mci>mi2_less));
        end
        sub=match(mci(mi2_less:mi2_more),:);
        md2=interp1(sub(:,4),sub(:,2),d2);
        newline1=[i c1 comp(a,3) md2 comp(a,5) d2 comp(a,7)];
        newline2=[i c2 md2 comp(a,4) d2 comp(a,6) new];
        
        % insert newline
        comp(a+2:end+1,:)=comp(a+1:end,:);
        comp(a,:)=newline1;
        comp(a+1,:)=newline2;
        section(ni+2:end+1,:)=[section(ni+1:end,1) section(ni+1:end,2:3)+1];
        newn=max(section(:,1))+1;
        section(ni+1,:)=[newn a+1 b+1];
        section(ni,:)=[n a a];
        %comp(a-1:b+2,:)
        %section(ni-1:ni+2,:)
        handles.comp=comp;
    else
        section(ni+2:end+1,:)=[section(ni+1:end,1) section(ni+1:end,2:3)];
        newn=max(section(:,1))+1;
        section(ni+1,:)=[newn a+1 b];
        section(ni,:)=[n a a];
        %disp('inserted')
        %section(ni-1:ni+2,:)
    end
    handles.section=section;
    guidata(gcbo, handles);
    axes(handles.axes1);
    x=xlim(handles.axes1);
    hold on;
    plot([0 x(2)],[d2 d2],'k:');
    text(handles.textx,d2,num2str(newn),'FontSize',7,'VerticalAlignment','top');
end

% --------------------------------------------------------------------
% Move section boundary or change its reference hole
function varargout = text_ButtondownFcn(h, eventdata, handles, varargin)

num=str2num(get(gcbo,'String'));
if(get(handles.edit,'Value')==1)
    targ=load(get(handles.sig1,'String'));
    hole{1}=targ;
    mhole{1}=hole{1};
    match{1}=hole{1};
    match{1}(:,3:4)=match{1}(:,1:2);
    maxi=2;
    for i=2:4
        sig=['get(handles.sig' num2str(i) ',''String'')'];
        if(~isempty(eval(sig)))
            maxi=i;
        else
            break;
        end
        msig=['get(handles.msig' num2str(i) ',''String'')'];
        m=['get(handles.match' num2str(i) ',''String'')'];
        hole{end+1}=load(eval(sig),'-ASCII');
        mhole{end+1}=load(eval(msig),'-ASCII');
        match{end+1}=load(eval(m),'-ASCII');
    end
    offset=handles.offset;  
    section=handles.section;
    comp=handles.comp;    
    ni=find(section(:,1)==num);
    a=section(ni,2);
    b=section(ni,3);
    oldi=comp(a,1);
    
    [x1,y1]=ginput(1);
    axes(handles.axes1)
    hold on;
    
    xx=xlim(handles.axes1);
    if(x1>.9*(xx(2)-xx(1))+xx(1))   % if clicked on right edge of plot, move segement boundary
        newi=oldi;
    else
        dist=[];
        % Estimate location of signals in plot to identify which hole was selected
        for i=1:maxi
            ind=find(mhole{i}(:,2)>=comp(a,5)-10 & mhole{i}(:,2)<=comp(b,6)+10);
            if(~isempty(ind))
                m=mean(mhole{i}(ind,3))+(i-1)*offset;
            else
                m=mean(mhole{i}(:,3))+(i-1)*offset;
            end
            dist(i)=abs(x1-m);
        end
        newi=find(dist==min(dist));
    end
    ind=[];
    
    %disp('START')
    %comp(section(ni,2):section(ni+1,3),:)
    
    if(newi==oldi)
        if(ni==1)
            disp('Error: Cannot move top of first segment.')
            return
        end
        a1=section(ni-1,2);
        b1=section(ni-1,3);
        newd=y1;
        if(newd<=comp(a1,5) | newd>=comp(b,6))
            disp('Error: Cannot cross section boundaries')
            comp(a1:b+3,:)
            section(ni-1:ni+2,:)
            return
        end
         x=get(gcbo,'Position');
         oldd=x(2);
         newx=x;
         newx(2)=newd;
%         set(gcbo,'Position',newx);          % move segment label 
%         delete(findobj('YData',[x(2) x(2)],'Parent',handles.axes1));
%         x2=xlim(handles.axes1);
%         plot([0 x2(2)],[newd newd],'k:');
        
        % if segment boundary moves up, change comp line that contains new boundary
        if(newd<comp(a,5))    
            h2=comp(a-1,1);   % hole number of previous segment, needs to change to h of current segment
            h=newi;           % hole number of current segment
            ind1=[a1:b1];     % comp indices for previous segment
            ii=find(comp(ind1,6)>=newd);  % all indices in comp that needs to change
            ind=ind1(ii);
            comp(ind(1),6)=newd;    % set end of previous segment to new boundary
            comp(ind(1)+1,5)=newd;  % comp(ind(1)+1,3)=depth in hole h will be set in for loop below
            c1=comp(ind(1),2);      % core number of previous segment
            %comp(ind(1)-1:ind(end)+1,:)
            
            % Find indices in match file for depth interpolation and transformation
            mci=find(match{h2}(:,1)>= c1-1 | match{h2}(:,1)<=c1+1);
            mi2_less=max(find(match{h2}(mci,4)<=newd & match{h2}(mci,3)==c1));
            mi2_more=min(find(match{h2}(mci,4)>newd & match{h2}(mci,3)==c1));
            if(isempty(mi2_less))
                mi2_less=max(find(match{h2}(mci,4)<=newd));
            end
            if(isempty(mi2_more) | mi2_more<=mi2_less)
                mi2_more=min(find(match{h2}(mci,4)>=newd & mci>mi2_less));
            end
            mi1_less=max(find(match{h}(:,4)<=newd));  % make sure mi1_less exists because will need it later
            if (isempty(mi1_less) | isempty(mi2_less) | isempty(mi2_more))
                disp('Error: Could not find requested depth in match file')
                return
            else  
                x=get(gcbo,'Position');
                oldd=x(2);
                newx=x;
                newx(2)=newd;
                set(gcbo,'Position',newx);          % move segment label 
                delete(findobj('YData',[x(2) x(2)],'Parent',handles.axes1));
                x2=xlim(handles.axes1);
                plot([0 x2(2)],[newd newd],'k:');
            end
            sub=match{h2}(mci(mi2_less:mi2_more),:);
            try
                md2=interp1(sub(:,4),sub(:,2),newd);
            catch
                newd
                md2=interp1(sub([1 end],4),sub([1 end],2),newd)
            end
            comp(ind(1),4)=md2;        % new segment end depth in previous hole 
            %disp('BOUNDARY UP. old and new boundary:')
            %[h2 oldd h newd]
            %comp(ind(1)-1:ind(1)+1,:)
            %section(ni-1:ni+1,:)
            
            section(ni-1,3)=ind(1);
            section(ni,2)=ind(1)+1;

            %oldsect=section(ni-1:ni+1,:);
            %section(ni-1:ni+1,:)
            ind2=ind;
            if(max(size(ind))>1)
                ind=ind(2:end);         % remaining comp indices that need to switch holes
            else
                ind=[];
            end
        % If segment boundary moves down, change comp line that contains new boundary
        elseif (newd>comp(a,5))     
            h=comp(a-1,1);          % hole number of previous segment
            h2=newi;                % hole number of current segment, needs to change to h of prev. segment
            ind1=[a:b];             % comp indices for current segment

            %comp(ind1,:)
            %newd
            ii=find(comp(ind1,5)<newd);     % all indices in comp that needs to change
            ind=ind1(ii);
            %comp(ind(1)-1:ind(end)+3,:)
            %section(ni-2:ni+1,:)
            if(ind(end)==b)
                comp(b+2:end+1,:)=comp(b+1:end,:);
                comp(b+1,:)=comp(b,:);
                section(ni+1:end,:)=[section(ni+1:end,1) section(ni+1:end,2:3)+1];
                section(ni,3)=section(ni,3)+1;
            end
                
            comp(ind(end),6)=newd;   % comp(ind(1)+1,4)=depth in hole h will be set in for loop below
            comp(ind(end)+1,5)=newd;
            c1=comp(ind(end),2);
            
            mci=find(match{h2}(:,1)>= c1-1 | match{h2}(:,1)<=c1+1);
            mi2_less=max(find(match{h2}(mci,4)<=newd & match{h2}(mci,3)==c1));
            mi2_more=min(find(match{h2}(mci,4)>newd & match{h2}(mci,3)==c1));
            if(isempty(mi2_less))
                mi2_less=max(find(match{h2}(mci,4)<=newd));
            end
            if(isempty(mi2_more) | mi2_more<=mi2_less)
                mi2_more=min(find(match{h2}(mci,4)>=newd & mci>mi2_less ));
            end
            mi1_more=max(find(match{h}(:,4)>=newd));  % make sure mi1_more exists because will need it later
            if (isempty(mi1_more) | isempty(mi2_less) | isempty(mi2_more))
                disp('Error: Could not find requested depth in match file')
                return
            else  
                x=get(gcbo,'Position');
                oldd=x(2);
                newx=x;
                newx(2)=newd;
                set(gcbo,'Position',newx);          % move segment label 
                delete(findobj('YData',[x(2) x(2)],'Parent',handles.axes1));
                x2=xlim(handles.axes1);
                plot([0 x2(2)],[newd newd],'k:');
            end

            sub=match{h2}(mci(mi2_less:mi2_more),:);
            try
                md=interp1(sub(:,4),sub(:,2),newd);
            catch
                newd
                md=interp1(sub([1 end],4),sub([1 end],2),newd)
            end
            
            %comp(ind(end),3)=md;      % new segment start depth in current hole 
            comp(ind(end)+1,3)=md;      % new segment start depth in current hole 

            %disp('BOUNDARY DOWN. old and new boundary:')
            %[h2 oldd h newd]
            %comp(ind(1)-1:ind(end)+1,:)
            %section(ni-1:ni+1,:)

            section(ni-1,3)=ind(end);
            section(ni,2)=ind(end)+1;
            
            %section(ni-1:ni+1,:)
            ind2=ind;
            ind=ind(1:end);
            ni=ni-1;
            %oldsect=section(ni-1:ni+1,:);
       end
    % If source segment changes, instead of a segment boundary move
    else            
        ind=[a:b];  
        ind2=ind;
        h=newi;
        h2=oldi;
        newd=comp(a,5);
        oldd=comp(b,6);
        
        % make sure depths found in match for new hole, otherwise return and display error
        mi1_less=max(find(match{h}(:,4)<=newd)); 
        mi2_less=max(find(match{h}(:,4)<=oldd));
        if (isempty(mi1_less) |  isempty(mi2_less))
            disp('Error: Could not find requested depth in hole')
            return
        end
        mind=[1:max(size(match{h}(:,1)))]';
        mi1_more=min(find(match{h}(mind,4)>=newd & mind>=mi1_less ));
        mi2_more=min(find(match{h}(mind,4)>=oldd & mind>=mi2_less ));
        if (isempty(mi1_more) | isempty(mi2_more))
            disp('Error: Could not find requested depth in hole')
            return
        end
    end
    %disp('handled boundary case')

    oldind=ind;
    
    if(~isempty(ind))      % change source hole in comp for either a boundary move or a source hole change
        %disp('comp before switching hole')
        %comp(ind,:)
        j=1;
        while j<=max(size(ind))
            %ind(j)
            d=comp(ind(j),5);           
            d2=comp(ind(j),6);
            k=min(find(targ(:,2)>=d));
            if isempty(k)
                disp('Error finding depth in target hole')
                return
            end
            c1=targ(k,1);        % first target core containing top of segment

            [md,mc1,c1a]=interpc(d,c1,[],match{h});
            %disp('Interpc results for d')
            %[c1a d mc1 md]
            if(isnan(md) | isnan(mc1) | isnan(c1a))
                disp('Error')
                d
                return
            end
            %[d md mc1] 
            
            [md2,mc2,c2]=interpc(d2,c1,mc1,match{h});
            %disp('Interpc results for d2')
            %[c2 d2 mc2 md2]
            if(isnan(md2) | isnan(mc2) | isnan(c2))
                disp('Error')
                d2
                return
            end
            %[d2 md2 mc2]
            
            if (mc1==mc2)
                oldline=comp(ind(j),:);
                newline=[h mc1 md md2 d d2];
                comp(ind(j),1:6)=[h mc1 md md2 d d2];
            else
                % if segment contains core gap add new entries to comp and section
                % assumes step contains only a single gap
                disp('Warning selected segment contains core gap')
                compi=ind(j);
                %comp(ind,:)
                %comp(compi,:)
                
                mtop=min(find(match{h}(:,1)==mc2));  % insert new entry in comp which starts at core top
                
                %if(match{h}(mtop,4)<=d2)
                
                mend=mtop-1;
                %match{h}(mend:mtop,:)
                
                % if no gap in target and segment cores overlap, switch core at end of first core
                if(match{h}(mtop,3)==match{h}(mend,3) & match{h}(mtop,4)<=match{h}(mend,4))
                    %disp('hole gap: case 1')
                    % find depth in mc2 which is matched to end of core mc1
                    [mtop2,mc2a,c2a]=interpc(match{h}(mend,4),c1a,mc2,match{h});
                    if(isnan(mtop2) | isnan(mc2a) | isnan(c2a))
                        disp('Error')
                        match(mend:mtop,:)
                        return
                    end
                    
                    % check core in target to set target-gap-flag
                    [mend1,mc1b,c1b]=interpc(match{h}(mend,4),c1a,match{h}(mend,1),match{h});
                    [md2b,mc2b,c2b]=interpc(d2,c1a,mc2,match{h});
                    if(isnan(c1b) | isnan(c2b))
                        disp('Error')
                        match{h}(mend:mtop,:)
                        return
                    end
                    flag1=(c1b>c1a);
                    flag2=(c2b>c1b);
                    
                    %comp: 1-hole, 2-core, 3-mbsf start, 4-mbsf stop, 5-target start, 6-target stop, 
                    %      7-flag for contains target coretop
                    oldline=comp(compi,:);
                    newline1=[h mc1   md match{h}(mend,2)  d  match{h}(mend,4)  flag1];
                    newline2=[h mc2a  mtop2           md2  match{h}(mend,4) d2  flag2];
                    
                else
                    %disp('hole gap: case 2')
                    if(comp(compi,7)==1)
                        if(match{h}(mtop,3)==c1)
                            flag1=0;
                            flag2=1;
                        else
                            flag1=1;
                            flag2=0;
                        end
                    else
                        flag1=0;
                        flag2=0;
                    end
                    
                    oldline=comp(compi,:);
                    newline1=[h mc1  md   match{h}(mtop,2)  d  match{h}(mtop,4)  flag1];
                    newline2=[h mc2  match{h}(mtop,2)  md2  match{h}(mtop,4) d2  flag2];
                end
                
                % insert newline
                comp(compi+2:end+1,:)=comp(compi+1:end,:);
                comp(compi,:)=newline1;
                comp(compi+1,:)=newline2;
                %increase j and the length of ind to account for insertion into comp
                %disp('update ind')
                %ind
                %compi
                j=j+1; 
                ind=[ind ind(end)+1];
                %ind
                %ind(j)

                %comp(section(ni-1,2):section(ni+1,3),:)
                %section(ni-1:ni+1,:)
                
                % section= segment label, start index in comp, end index in comp
                section(ni+1:end,:)=[section(ni+1:end,1) section(ni+1:end,2:3)+1];
                section(ni,3)=section(ni,3)+1;

                %section(ni-1:ni+1,:)
                %compi=compi+1;
            end
            j=j+1;    
        end   %end while loop
        
        %oldind
        %ind
        %oldsect
        %section(ni-1:ni+1,:)
        
    end       %end if(~isempty(ind))
    
    %disp('END:')
    %comp(section(ni-1,2):section(ni+1,3),:)
    %if(~isempty(ind))
    %    comp(ind(1)-1:ind(end)+2,:)
    %end
    
    handles.comp=comp;
    handles.section=section;
    guidata(gcbo, handles);
   
    d1=min([newd oldd]);
    d2=max([newd oldd]);
    ind=find(mhole{h2}(:,2)>=d1 & mhole{h2}(:,2)<=d2);
    color=handles.color;
    plot(mhole{h2}(ind,3)+(h2-1)*handles.offset,mhole{h2}(ind,2),'Color',color(h2,:));
    ind=find(mhole{h}(:,2)>=d1 & mhole{h}(:,2)<=d2);
    plot(mhole{h}(ind,3)+(h-1)*handles.offset,mhole{h}(ind,2),'r');
    
end



% --------------------------------------------------------------------
% Interpolate from depth in target to depth in hole h
function [mx, mc1, c1] = interpc(x, c, mc, match)

%[x c mc]
mc_in=~isempty(mc);
err=0;
if isempty(mc)
    try
        mc=match(min(find(match(:,4)>=x & match(:,3)==c)),1);
    catch
        err=1;
    end
    if(isempty(mc) | err)
        ind=min(find(match(:,4)>=x));
        if ~isempty(ind)
            mc=match(ind,1);
            c=match(ind,3);
        else
            s=['Error: Depth ' x ' not found in target'];
            disp(s)
            mx=NaN;
            mc1=NaN;
            c1=NaN;
            return
        end
    end
end

sub=find(match(:,1)==mc & match(:,3)==c);
%x
%match(sub,:)

if (isempty(sub) | min(match(sub,4)>x))
    f=max(find(match(:,4)<=x));
    sub=[f; sub];
end
if (max(match(sub,4)<x))
    g=min(find(match(:,4)>=x));
    jnd=find(match(:,1)==match(g,1) & match(:,3)==match(g,3));
    if(min(match(jnd,4))<=x)
        sub=jnd;
    else
        sub=[sub; g];
    end
end

%disp('sub after limit check')
%match(sub,:)

i=find(match(sub,4)==x);
if(~isempty(i))
    mx=match(sub(i(1)),2);
else
    try 
        mx=interp1(match(sub,4),match(sub,2),x);
    catch
        try
            mx=interp1(match(sub([1 end]),4),match(sub([1 end]),2),x);
        catch
            s=['Error: Depth ' num2str(x) ' not found in target'];
            disp(s)
            mx=NaN;
            mc1=NaN;
            c1=NaN;
            return
        end
    end
end
i=max(find(match(sub,4)<=x));
%[match(sub(i),:) x]

% if mx is at a core top and the calling function requested the previous core, 
% return the previous core info so it can be used as a segment end in comp
if(match(sub(i),4)==x & sub(i)>1 & match(sub(i)-1,1)==mc & mc_in==1)
    %disp('interpc special case')
    %match(sub(i)-1:sub(i),:)
    sub(i)=sub(i)-1;
end
mc1=match(sub(i),1);
c1=match(sub(i),3);
return


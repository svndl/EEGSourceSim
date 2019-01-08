function h = textfit(x,y,txt,varargin)
% textfit(x,y,txt,varargin)
%
% Mike Lawrence 2011

ythresh = 0.4;  % maximal allowable overlap (includes cell padding accounted for in "extent" property)
xthresh = 0.1;

n = length(x);
if ~iscell(txt), txt={txt}; end
if length(y)~=n || length(txt)~=n, error('length mismatch between x,y,txt'); end

h = text(x,y,txt,varargin{:});

yl=ylim; ytot=diff(yl);
xl=xlim; xtot=diff(xl);

maxtries = 100;
for t=1:maxtries
  ext = nan(n,4);
  for i=1:n, ext(i,:) = get(h(i),'extent'); end
  xstart=ext(:,1); xsz=ext(:,3); xend=xstart+xsz;
  ystart=ext(:,2); ysz=ext(:,4); yend=ystart+ysz;
  overlapx = zeros(n,n);
  overlapy = zeros(n,n);
  for i1=1:n-1, for i2=i1+1:n
    if xstart(i1)<=xend(i2)&xstart(i2)<=xend(i1)
      overlapx(i1,i2)=(min(xend(i2)-xstart(i1),xend(i1)-xstart(i2)))/(min(xsz(i1),xsz(i2)));
    end
    if ystart(i1)<=yend(i2)&ystart(i2)<=yend(i1)
      overlapy(i1,i2)=(min(yend(i2)-ystart(i1),yend(i1)-ystart(i2)))/(min(ysz(i1),ysz(i2)));
    end
  end,end
  overlapmax = max(overlapx,overlapy);
  ov = (overlapx>xthresh & overlapy>ythresh);
  [o1 o2] = find(ov);
  if isempty(o1), break; end
  [tmp ord] = sort(overlapmax(find(ov)));
  o1=o1(ord); o2=o2(ord);
  moved = false(n,1);
  for i=1:length(o1), i1=o1(i); i2=o2(i);
    if moved(i1) || moved(i2), continue; end
    pos1 = get(h(i1),'position');
    pos2 = get(h(i2),'position');
    oy = overlapy(i1,i2)*min(ysz(i1),ysz(i2));
    ox = overlapx(i1,i2)*min(xsz(i1),xsz(i2));
    if oy/ytot < ox/xtot   % overlapy is easier to fix
      shift = 0.5*(1-ythresh)*oy;
      if ystart(i1)<ystart(i2)  % i1 above i2
        pos1(2)=pos1(2)-shift; pos2(2)=pos2(2)+shift;
      else                      % i1 below i2
        pos1(2)=pos1(2)+shift; pos2(2)=pos2(2)-shift;
      end
    else                                 % overlapx is easier to fix
      shift = 0.5*(1-xthresh)*ox;
      if xstart(i1)<xstart(i2)  % i1 left of i2
        pos1(1)=pos1(1)-shift; pos2(1)=pos2(1)+shift;
      else                      % i1 right of i2
        pos1(1)=pos1(1)+shift; pos2(1)=pos2(1)-shift;
      end
    end
    set(h(i1),'position',pos1);
    set(h(i2),'position',pos2);
    moved([i1 i2]) = true;
  end
end

if nargout==0, clear h, end
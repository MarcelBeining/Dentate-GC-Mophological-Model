% FUNCTION interp_border interpolates the contour "line" either by
% distributing a defined amount of points in equidistance (mode 1),in a defined
% distance (mode 2) or dependent on a point distribution histogram. Additionally, smoothing can be applied to the contour.
% if no pdistr is given, a gaussian distribution is shown

function outline = interp_border(line,x,mode,options,filtsize,pdistr,min_intervall)

if nargin < 4 || isempty(options)
    options = '';
end
if mode >= 3 && (nargin < 6 || isempty(pdistr))
    pdistr = -100:100;
    pdistr = 10/(20*sqrt(2*pi)) * exp(-0.5*(pdistr/(2*50)).^2);
end
if mode >= 3 && (nargin < 7 || isempty(min_intervall))
    min_intervall = Inf;
end

outline = [];
if size(line,1) < 2  %too low number of points
    return
end


d_line = diff(line,1,1);    % calc vectors from one to next point
lengths = sqrt(sum(d_line.^2,2));   % calc lengths of these vectors
d_line = d_line./repmat(lengths,[1 size(d_line,2)]);    % norm the vectors
glength = sum(lengths); % calc length of whole contour
if mode == 1
    if x < 2    %too low number of points
        return 
    end
    pnumber = x;
    ilength = glength/(pnumber-1);  % calc distance between interpolated points
elseif mode == 2
    if x > glength  % too large intervall
        return
    end
    ilength = x;
    pnumber = round(glength/ilength)+1; % calc number of intervall points 
    % CAUTION: End point of interpolated line may change drastically if
    % intervall length is huge!
elseif mode >= 3
    pnumber = x;
    pmindistr = ones(1,floor(glength)) * (1/(min_intervall));%*(pnumber-2)));%./ (pnumber-2);
    if sum(pmindistr) > pnumber-1
        pmindistr = pmindistr * 0;      %min_intervall will be ignored if too many points would be set
    end
%     if numel(pdistr) ~= pnumber-1     %interpolate distribution if it does not fit to number of points
        rat = (numel(pdistr)-1)/(floor(glength)-1);
        pdistr = interp1(1:numel(pdistr),pdistr,1:rat:numel(pdistr),'pchip');
        if mode == 4
            pdistr(1) = 0;
        end
%     end
%     pdistr = -pdistr + max(pdistr)+1; % invert the distribution so that huge values will give narrow points..the +1 is necessary that not two points are at the same position
       pdistr = pdistr / sum(pdistr)*((pnumber-1)-sum(pmindistr)) + pmindistr;
%     pdistr = pdistr / sum(pdistr);  % norm it again
else
    return
end
if mode >= 3
%     cum_ilength = cumsum(pdistr*glength);
%     cumsetpoint = cumsum(pdistr*(pnumber-1));

    cum_ilength = interp1(cumsum(pdistr),1:numel(pdistr),1:pnumber-1);
else
    cum_ilength = cumsum(repmat(ilength,[pnumber-1,1]));    %calc cumulative sum of these distances
end

cum_ilength(end) = glength;     % avoid calc-error when cumulating ilengths

cum_lengths = cumsum([0;lengths]);      %calc distances of original points to first point (along contour)

outline = zeros(pnumber,size(line,2));  % predefine output array
outline(1,:) = line(1,:);   %first point is original point
% 
% if mode == 3
%     
%     c=1;
%     p=1;
%     while p < pnumber-1
%         if cumsetpoint(c) < p
%             c=c+1;
%         else
%             pdiff = cumsetpoint(c);
%             if p == pnumber -2
%                 'g'
%             end
%             while pdiff > p 
%                 %% THERE IS AN END PROBLEM =(
%                 l = interp1(cumsetpoint([c-1,c]),[c-1 c],p);
%                 ind = find(cum_lengths < l,1,'last');
%                 outline(p,:) = line(ind,:)+(line(ind+1,:)-line(ind,:))/sqrt(sum((line(ind+1,:)-line(ind,:)).^2))*(l-cum_lengths(ind));
%                 pdiff = pdiff - 1;
%                 p=p+1;
%             end
%         end
%     end
% else
    for p = 1:pnumber-1
        ind = find(cum_ilength(p)>cum_lengths,1,'last');    % find the original point which is closest to the current distance
        outline(p+1,:) = line(ind,:) + d_line(ind,:)*(cum_ilength(p)-cum_lengths(ind));   %calc position of interpolated point by using the closest original point and the vector to next point
    end
% end

if strfind(options,'smooth')
    %     filtsize = 0.05;
    if round(pnumber*filtsize) ~= pnumber*filtsize
        kern = rem(floor(pnumber*filtsize),2) * floor(pnumber*filtsize) + rem(ceil(pnumber*filtsize),2) * ceil(pnumber*filtsize); % rounds to the odd number
    else
        kern = pnumber*filtsize + ~rem(pnumber*filtsize ,2) ;
    end
    tmp = outline([1,end],:);
    outline(:,1) = convn (padarray(outline(:,1),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid');
    outline([1,end],:) = tmp;   % accounts for smoothing error at edges
    outline(:,1) = flipud( convn (padarray(flipud(outline(:,1)),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid'));
    outline(:,2) = convn (padarray(outline(:,2),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid');
    outline([1,end],:) = tmp;   % accounts for smoothing error at edges
    outline(:,2) = flipud(convn (padarray(flipud(outline(:,2)),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid'));
    outline([1,end],:) = tmp;   % accounts for smoothing error at edges
end
function pts_GCL = giveGCLpoints3(N_GCL,GCLdens,Somata,thicknesses,pts_IML,shapeprofile)
show = 0;
fac = 2;
while 1
    if sum(pts_IML(:,3) < sum(thicknesses(1:2))+thicknesses(3)/fac) == 0
        fac = fac/2;
    else
        break
    end
end
pts_IML = pts_IML(pts_IML(:,3) < sum(thicknesses(1:2))+thicknesses(3)/fac,:);  % only use IML points in the first half of IML as reference to GCL
refpoint = sum(thicknesses(1:2));
maxdist = 0;
% if show
    pts_IMLorig = pts_IML;
% end
for n = 1:size(pts_IML,1) % scale all ML points to reference level
    pts_IML(n,1) = interp1([Somata(3),pts_IML(n,3)],[Somata(1),pts_IML(n,1)],refpoint);
    pts_IML(n,2) = interp1([Somata(3),pts_IML(n,3)],[Somata(2),pts_IML(n,2)],refpoint);
    maxdist = max(maxdist,sqrt(sum((pts_IML(n,1:2)-Somata(1:2)).^2))); % find point that is most far..approximates the cone radius
end
nclusters = 1;
while 1
    flag = false;
    [ind,centroids,~,dist] = kmeans(pts_IML(:,1:2),nclusters); % note: extracting distance from this function returns quadratic distances!
    %         sumD ./ hist(ind,1:max(ind))
    maxclustdist = zeros(nclusters,1);
    for n = 1:nclusters
        %         dist = sqrt(sum((pts_IML(ind==nclusters,1:2) - repmat(centroids(n,:),sum(ind==nclusters),1)).^2,2));  % calculate all distances from centroid
        maxclustdist(n) = max(sqrt(dist(ind==n,n)));
        if  maxclustdist(n) > maxdist/2 && nclusters < 2  % if inner cluster distance is greater than half the cone radius, increase amount of clusters
            flag = true;
            break
        end
    end
    if flag
        nclusters = nclusters + 1; % increase amount of kmeans clusters
    else
        break
    end
end
% display(nclusters)
pointcount = hist(ind,1:nclusters);
pointcount = cumsum(pointcount/sum(pointcount));
if show
    figure;hold all,plot3(Somata(1),Somata(2),Somata(3),'Marker','^','Color','r','LineWidth',5)
    col = colorme(nclusters);
    for n = 1:nclusters
        plot3(centroids(n,1),centroids(n,2),refpoint,'Marker','o','Color',col{n},'LineWidth',2)
        %     plot3(pts_IML(ind==n,1),pts_IML(ind==n,2),repmat(refpoint,sum(ind==n),1),'Marker','x','Color',col{n},'LineStyle','None')
        plot3(pts_IMLorig(ind==n,1),pts_IMLorig(ind==n,2),pts_IMLorig(ind==n,3),'Marker','x','Color',col{n},'LineStyle','None')
    end
end
% N_GCL = max(N_GCL,nclusters); % have at least as many gcl points as clusters
if N_GCL > 0
    [~, z] = find(diff(cat(2,zeros(N_GCL,1),repmat(rand(N_GCL,1),1,numel(GCLdens)) <= repmat(GCLdens,N_GCL,1)),1,2));
    zindices = Somata(3):(sum(thicknesses(1:2))-Somata(3))/(numel(GCLdens)-1):sum(thicknesses(1:2));  % scale IML thickness to GCLdens (normally 100 points)
    z = zindices(z);
    
%     scale = normrnd(0,0.25,1,N_GCL); % calculate the random scaling of the centroid radius
%     scale(abs(scale)>1) = 1;
%     x = scale .* sin(rand(1,N_GCL)*2*pi);  % calculate the relative x positions somewhere in a circle
%     y = scale .* cos(rand(1,N_GCL)*2*pi); % calculate the relative x positions somewhere in a circle
    x = normrnd(0,0.2,1,N_GCL); % calculate the random scaling of the centroid radius
    y = normrnd(0,0.2,1,N_GCL); % calculate the random scaling of the centroid radius
    xorig = x;
    yorig = y;
    len = sqrt(x.^2+y.^2);
    x(len>1) = x(len>1)/len(len>1);
    y(len>1) = y(len>1)/len(len>1);
    
    thisconex = interp1(shapeprofile(:,3)+Somata(3),shapeprofile(:,1),cat(2,z,refpoint),'linear','extrap');
    thisconey = interp1(shapeprofile(:,3)+Somata(3),shapeprofile(:,2),cat(2,z,refpoint),'linear','extrap');
    
    for n = 1:N_GCL
        if n <= nclusters
            tind = n; % first give every cluster at one IML point
        else
            %             tind = randi(nclusters);  % random cluster choose
            tind = find(rand<=pointcount,1,'first'); % clusters with more points in it are more likely to get more IML points
        end
        newcentroid(1) = (centroids(tind,1)-Somata(1)) / thisconex(end) * thisconex(n) + Somata(1);
        newcentroid(2) = (centroids(tind,2)-Somata(2)) / thisconey(end) * thisconey(n) + Somata(2);
        newmaxclustdist = maxclustdist(tind) / max([thisconex(end)/thisconex(n),thisconey(end)/thisconey(n)]);  % this is an approximation because I am too lazy
        
        x(n) = x(n) * newmaxclustdist + newcentroid(1);
        y(n) = y(n) * newmaxclustdist + newcentroid(2);
        if show
            plot3(x(n),y(n),z(n),'Marker','d','Color',col{tind},'LineWidth',2)
        end
    end
    
    pts_GCL = [x',y',z'];
    
else
    pts_GCL = zeros(0,3);
end
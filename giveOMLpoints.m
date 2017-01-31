function pts_OML = giveOMLpoints(N_OML,OMLdens,Somata,thicknesses,pts_ML,shapeprofile)
show = 0;
% fac = 2;
% while 1
%     if sum(pts_ML(:,3) < sum(thicknesses(1:3))+sum(thicknesses(4:5))/fac) == 0
%         fac = fac/2;
%     else
%         break
%     end
% end
% pts_ML = pts_ML(pts_ML(:,3) < sum(thicknesses(1:3))+sum(thicknesses(4:5))/fac,:);  % only use IML points in the first half of IML as reference to GCL

refpoint = sum(thicknesses(:));
maxdist = 0;
% if show
    pts_MLorig = pts_ML;
% end
for n = 1:size(pts_ML,1) % scale all ML points to reference level
    pts_ML(n,1) = interp1([Somata(3),pts_ML(n,3)],[Somata(1),pts_ML(n,1)],refpoint,'linear','extrap');
    pts_ML(n,2) = interp1([Somata(3),pts_ML(n,3)],[Somata(2),pts_ML(n,2)],refpoint,'linear','extrap');
    maxdist = max(maxdist,sqrt(sum((pts_ML(n,1:2)-Somata(1:2)).^2))); % find point that is most far..approximates the cone radius
end
nclusters = 1;
while 1
    flag = false;
    [ind,centroids,~,dist] = kmeans(pts_ML(:,1:2),nclusters); % note: extracting distance from this function returns quadratic distances!
    %         sumD ./ hist(ind,1:max(ind))
    maxclustdist = zeros(nclusters,1);
    for n = 1:nclusters
        %         dist = sqrt(sum((pts_ML(ind==nclusters,1:2) - repmat(centroids(n,:),sum(ind==nclusters),1)).^2,2));  % calculate all distances from centroid
        maxclustdist(n) = max(sqrt(dist(ind==n,n)));
        if  maxclustdist(n) > maxdist/1.5  % if inner cluster distance is greater than half the cone radius, increase amount of clusters
            flag = true;
            break
        end
    end
%     if nclusters == N_OML
%         break
%     end
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
        %     plot3(pts_ML(ind==n,1),pts_ML(ind==n,2),repmat(refpoint,sum(ind==n),1),'Marker','x','Color',col{n},'LineStyle','None')
        plot3(pts_MLorig(ind==n,1),pts_MLorig(ind==n,2),pts_MLorig(ind==n,3),'Marker','x','Color',col{n},'LineStyle','None')
    end
end
% N_OML = max(N_OML,nclusters); % have at least as many iml points as clusters
if N_OML > 0
    [~, z] = find(diff(cat(2,zeros(N_OML,1),repmat(rand(N_OML,1),1,numel(OMLdens)) <= repmat(OMLdens,N_OML,1)),1,2));
    zindices = sum(thicknesses(1:4)):thicknesses(5)/(numel(OMLdens)-1):sum(thicknesses(1:5));  % scale IML thickness to IMLdens (normally 100 points)
    z = zindices(z);
%     sum((z-sum(thicknesses(1:4)))< 7)
    z((z-sum(thicknesses(1:4)))< 7) = sum(thicknesses(1:4))+7;  % 7µm Sperrzone sonst rutschen BPs in den Layer darunter mit smooth tree..
    scale = normrnd(0,0.1,1,N_OML); % calculate the random scaling of the centroid radius
    scale(abs(scale)>1) = 1;
    p = rand(1,N_OML);
    x = scale .* sin(p*2*pi);  % calculate the relative x positions somewhere in a circle
    y = scale .* cos(p*2*pi); % calculate the relative x positions somewhere in a circle
    
    thisconex = interp1(shapeprofile(:,3)+Somata(3),shapeprofile(:,1),cat(2,z,refpoint),'linear','extrap');
    thisconey = interp1(shapeprofile(:,3)+Somata(3),shapeprofile(:,2),cat(2,z,refpoint),'linear','extrap');
    
    for n = 1:N_OML
        if n <= nclusters
            tind = n; % first give every cluster at one IML point
        else
            %             tind = randi(nclusters);  % random cluster choose
            tind = find(rand<=pointcount,1,'first'); % clusters with more points in it are more likely to get more IML points
        end
        newcentroid(1) = (centroids(tind,1)-Somata(1)) / thisconex(end) * thisconex(n) + Somata(1);
        newcentroid(2) = (centroids(tind,2)-Somata(2)) / thisconey(end) * thisconey(n) + Somata(2);
        newmaxclustdist = maxclustdist(tind) / min([thisconex(n)/thisconex(end),thisconey(n)/thisconey(end)]);  % this is an approximation because I am too lazy
        
        x(n) = x(n) * newmaxclustdist + newcentroid(1);
        y(n) = y(n) * newmaxclustdist + newcentroid(2);
        if show
            plot3(x(n),y(n),z(n),'Marker','d','Color',col{tind},'LineWidth',2)
        end
    end
    
    pts_OML = [x',y',z'];
    
else
    pts_OML = zeros(0,3);
end
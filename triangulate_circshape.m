function faces = triangulate_circshape(nu,nv)
% nv: number of circular ordered points per level
% nu: number of height levels

A = repmat(cat(2,(1:nv)',(1:nv)'+nv,circshift((1:nv)',-1)),nu,1); % first clockwise triangles of 4points
B = repmat(cat(2,(1:nv)',circshift((1:nv)',-1),circshift((1:nv)',-1)-nv),nu,1); %second clockwise triangle of 4 points
add = repmat(reshape(repmat(0:nu-1,nv,1),(nu)*(nv),1),1,3) * nv; % the different levels of the cone for first triangle set
A = A + add;  % add it
add = repmat(reshape(repmat(1:nu,nv,1),(nu)*(nv),1),1,3) * nv;% the different levels of the cone for second triangle set
B = B + add;  % add it
faces = cat(1,A,B);
faces(faces(:,2) > nv*nu,:) = []; % delete non existing triangles
faces(faces(:,3) < 1,:) = [];   % delete non existing triangles
faces = cat(1,faces,cat(2,ones(nv-2,1),(2:nv-1)',(3:nv)'),cat(2,ones(nv-2,1),(3:nv)',(2:nv-1)')+(nv*(nu-1)));  % add the deckels
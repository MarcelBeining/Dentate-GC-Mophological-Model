% XMST_TREE   Extended Minimum spanning tree based tree constructor.
% (trees package)
%
% [tree indx] = xxMST_tree (msttrees, X, Y, Z, L, nsyn, bf, thr, DIST,
% options)
% -------------------------------------------------------------------------
%
% Creates trees corresponding to the minimum spanning tree keeping the path
% length to the root small (with balancing factor bf). A sparse distance
% matrix DIST between nodes is added to the cost function. Don't forget to
% include input tree nodes into the distance matrix DIST. Input nodes may
% be grouped to labels for which if one node from that label is used all
% others from the same label are removed.
% For speed and memory considerations an area of close vicinity is drawn
% around each tree as it grows.
%
% Input
% -----
% - msttrees::vector: indices to the starting points of trees(# determines
%     # of trees), or starting trees as cell array of trees structures
%     {DEFAULT: additional node (0,0,0)}
% - X::vertical vector: X coords of pts to be connected {1000 rand. pts}
% - Y::vertical vector: Y coords of pts to be connected {1000 rand. pts}
% - Z::vertical vector: Z coords of pts to be connected {DEFAULT: zeros}
% - L::vertical vector: index attributing each node to a different label
% - bf::number between 0 1: balancing factor {DEFAULT: 0.4}
% - thr::value: max distance that a connection can span {DEFAULT: 20}
% - DIST::sparse matrix BIGNxBIGN: zero indicates probably no connection,
%     numbers increasing probabilities of a connection {DEFAULT: sparse
%     zeros matrix}. order of elements is first all trees in order and then
%     all open points.
% - maxT::number: maximum number of stems from root
% - options::string: {DEFAULT '-w'}
%     '-s' : show plot (much much much slower)
%     '-w' : with waitbar
%     '-t' : time lapse save
%     '-b' : suppress multifurcations
%     '-p' : relative paths
%     '-i' : add points in a tree-by-tree fashion
%
% Output
% ------
% if no output is declared the trees are added in trees
% - tree:: structured output trees, cell array if many
% - indx:: index indicating where points ended up [itree inode]
%
% See also rpoints_tree quaddiameter_tree BCT_tree MST_tree
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function [tree, indx] = xxMST_tree (msttrees, X, Y, Z, L, nsyn, bf, ...
    thr, DIST, maxT, plthr, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1) || isempty (msttrees)
    % starting tree is just point (0, 0, 0)
    msttrees       =           {};
    msttrees{1}.X  =            0;
    msttrees{1}.Y  =            0;
    msttrees{1}.Z  =            0;
    msttrees{1}.dA =   sparse (0);
    msttrees{1}.D  =            1;
    msttrees{1}.R  =            1;
    msttrees{1}.rnames = {'tree'};
end

if (nargin < 2) || isempty (X)
    X        = rand  (1000, 1)        .* 400;
end

if (nargin < 3) || isempty (Y)
    Y        = rand  (size (X, 1), 1) .* 400;
end

if (nargin < 4) || isempty (Z)
    Z        = zeros (size (X, 1), 1);
end

if (nargin < 5) || isempty (L)
    L        = [];
end

if (nargin < 6) || isempty (nsyn)
    nsyn     = 1;
end

if (nargin < 7) || isempty (bf)
    bf       = 0.4;
end

if (nargin < 8) || isempty (thr)
    thr      =  20;
end

if (nargin < 9) || isempty (DIST)
    DIST     =  [];
else
    Derr     = max (max (DIST));
    DIST     = DIST / Derr;
end

if (nargin < 10) || isempty (maxT)
    maxT     =  [];
end
if (nargin < 11) || isempty (plthr)
    plthr     =  [];
end
if (nargin < 12) || isempty (options)
    options  = '-w';
end

if ~iscell (msttrees)
    ID       = msttrees;
    msttrees = cell (1, length (ID));
    for ward = 1 : length (ID)
        msttrees {ward}.X      = X (ID (ward));
        msttrees {ward}.Y      = Y (ID (ward));
        msttrees {ward}.Z      = Z (ID (ward));
        msttrees {ward}.dA     =    sparse (0);
        msttrees {ward}.D      =             1;
        msttrees {ward}.R      =             1;
        msttrees {ward}.rnames =      {'tree'};
    end
    X (ID)   = [];
    Y (ID)   = [];
    Z (ID)   = [];
    if ~isempty (L)
        L (ID) = [];
    end
end

lenX         = length (X);        % number of targets
lent         = length (msttrees); % number of trees

if ~isempty  (DIST)
    iDIST    = cell (lent, 1);
    iDISTP   = cell (lent, 1);
    TSUM     = 0;
    for ward = 1 : lent,
        N    = length (msttrees {ward}.X); % number of nodes in tree;
        % DIST index creation, which indicates which node in the tree
        % corresponds to which field in DIST:
        iDIST {ward} = TSUM + 1 : TSUM + N;
        TSUM         = TSUM + N;
    end
end

if strfind (options, '-t')       % time lapse save
    timetrees    = cell (lent, 1);
    for tissimo  = 1 : lent
        timetrees{tissimo}{1} = msttrees {tissimo};
    end
end

if strfind (options, '-s')       % prepare a plot if showing
    % choose colors for the different trees:
    colors = [   ...
        [1 0 0]; ... 
        [0 1 0]; ... 
        [0 0 1]; ... 
        [0.2 0.2 0.2]; ... 
        [1 0 1]; ... 
        [1 1 0]; ... 
        [0 1 1]];
    if lent > 7
        colors = [colors; (rand (lent - 7, 3))];
    end
    clf; hold on;
    plot3    (X, Y, Z, 'k.');    % plot all points
    HP       = cell (1, lent);   % plot all starting trees
    for ward = 1 : lent
        HP {ward} = plot_tree (msttrees {ward});
    end
    view     (2);
    grid     on;
    axis     image;
end

% initialization:
N            =  cell (lent, 1);  % number of nodes in each tree
tthr         =  cell (lent, 1);  % vicinity threshold for each tree
root_dist    =  cell (lent, 1);  % dist. from targets to root in each tree
rdist        =  cell (lent, 1);
irdist       =  cell (lent, 1);
plen         =  cell (lent, 1);  % path length to the root within each tree
plencost     =  cell (lent, 1);  % path length cost relative to euclidean
avic         =  cell (lent, 1);
inX          =  cell (lent, 1);
dplen        =  cell (lent, 1);
ITREE        =  cell (lent, 1);
for ward     = 1 : lent
    N {ward} = length (msttrees {ward}.X); % number of nodes in tree
    if N {ward} > 1  % starting tree is not empty
        % threshold distance determining the vicinity circle:
        eucl        = eucl_tree (msttrees {ward});
        tthr {ward} = max (eucl) + thr;
        % calculate distance from all targets to root
        root_dist {ward} = sqrt ( ...
            (X - msttrees{ward}.X(1)).^2 + ...
            (Y - msttrees{ward}.Y(1)).^2 + ...
            (Z - msttrees{ward}.Z(1)).^2)';
        % starting path length to the root in the tree:
        plen {ward} = Pvec_tree (msttrees{ward});
        if strfind (options, '-p')  % relative path length cost
            plencost {ward} = plen{ward} - root_dist{ward};
        else
            plencost {ward} = plen {ward};
        end
        % calculate distance from all targets to any point on the tree
        dis  = zeros (1, lenX);
        idis = ones  (1, lenX);
        if strfind (options, '-b') % avoid multifurcations
            % find non-branch points:
            iCT = find (sum (msttrees{ward}.dA, 1) < 2);
            if ~isempty (maxT) % there is a max number of stems at root
                if sum (msttrees{ward}.dA (:, 1), 1) > maxT - 1
                    iCT (iCT == 1) = [];
                else
                    iCT = unique ([1 iCT]);
                end
            end
            % search only among non-branch-points:
            for te   = 1 : lenX
                sdis = sqrt ( ...
                    (X (te) - msttrees {ward}.X (iCT)).^2 + ...
                    (Y (te) - msttrees {ward}.Y (iCT)).^2 + ...
                    (Z (te) - msttrees {ward}.Z (iCT)).^2);
                % dis contains closest non-branch node on tree:
                [dis(te) idis(te)] = min (sdis);
            end
            idis     = iCT (idis); % retranslate index to all nodes
        else
            for te   = 1 : lenX
                sdis = sqrt ( ...
                    (X (te) - msttrees {ward}.X).^2 + ...
                    (Y (te) - msttrees {ward}.Y).^2 + ...
                    (Z (te) - msttrees {ward}.Z).^2);
                % dis contains closest node on tree:
                [dis(te) idis(te)] = min (sdis);
            end
        end
        % this was actually wrong because then points that are too far away
        % from current tree never gonna be considered!!!
%         % don't allow to go beyond the threshold distance:
%         dis (dis > thr) = NaN;
%         dis (dis > plthr) = NaN;
        % sort points according to their distance to the tree:
        [rdist{ward} irdist{ward}] = sort (dis);
        % set actual vicinity to all points in distance tthr of root:
        avic{ward} = sum (rdist{ward} < tthr{ward});
        if strfind (options, '-s'),
            plot3 ( ...
                X (irdist {ward} (1 : avic {ward})), ...
                Y (irdist {ward} (1 : avic {ward})), ...
                Z (irdist {ward} (1 : avic {ward})), 'g.');
        end
        % vector index in XYZ all points which are in vicinity but not yet
        % on tree:
        inX {ward} = irdist {ward} (1 : avic {ward});
        if ~isempty (DIST),
            % index of open points in distance matrix DIST:
            iDISTP {ward} = inX {ward} + TSUM;
            % initialize distance vector including path to root and extra
            % distance:
            dplen {ward}  = ...
                rdist {ward} (1 : avic {ward}) + ...
                bf * plencost {ward} (idis (inX {ward}))' + ...
                Derr * (1 - DIST ( ... % extra error term from DIST matrix
                iDIST  {ward} (idis (inX {ward}))', ...
                iDISTP {ward}));
        else
            % initialize distance vector including path to root
            dplen {ward}  = ...
                rdist {ward} (1 : avic {ward}) + ...
                bf * plencost {ward} (idis (inX {ward}))';
        end
        % initialize index vector indicating to which point on tree an open
        % point is closest to:
        ITREE {ward} = idis (inX {ward});
    else           % if tree is empty (this is easier)
        % starting path length to the root in the tree is just 0
        plen     {ward} =   0;
        plencost {ward} =   0;
        % threshold distance determining the vicinity circle
        tthr     {ward} = thr;
        % calculate distance from all open points to root
        root_dist{ward} = sqrt ( ...
            (X - msttrees{ward}.X(1)).^2 + ...
            (Y - msttrees{ward}.Y(1)).^2 + ...
            (Z - msttrees{ward}.Z(1)).^2)';
        % dis contains closest node on tree
        % this is simply the distance to root in this case
        dis  = root_dist{ward};
        % this was actually wrong because then points that are too far away
        % from current tree never gonna be considered even when tree grows there!!!
%         % don't allow to go beyond the threshold distance:
%         dis (dis > thr) = NaN;
%         dis (dis > plthr) = NaN;
        % sort points according to their distance to the root:
        [rdist{ward} irdist{ward}] = sort(root_dist{ward});
        % set actual vicinity to all points in distance tthr of root:
        avic{ward} = sum (rdist{ward} < tthr{ward});
        if strfind (options, '-s')
            plot3 ( ...
                X (irdist {ward} (1 : avic {ward})), ...
                Y (irdist {ward} (1 : avic {ward})), ...
                Z (irdist {ward} (1 : avic {ward})), 'g.');
        end
        % vector index in XYZ all points which are in vicinity but not yet
        % on tree:
        inX{ward} = irdist {ward} (1 : avic {ward});
        if ~isempty (DIST),
            % index of open points in distance matrix DIST:
            iDISTP {ward} = inX {ward} + TSUM;
            % initialize distance vector including path to root and extra
            % distance from DIST matrix:
            dplen{ward} = dis (inX {ward}) + ...
                Derr * (1 - DIST (iDIST{ward}(1),iDISTP{ward}));
        else
            % initialize distance vector including path to root
            dplen{ward} = dis(inX {ward});
        end
        % initialize index vector indicating to which point on tree an open
        % point is closest to. In the beginning all points are closest to
        % the root (#1):
        ITREE{ward} = ones (1, avic {ward});
    end
end

if strfind   (options, '-w')
    HW       = waitbar (0, 'finding minimum spanning tree...');
    set      (HW, 'Name', '..PLEASE..WAIT..YEAH..');
end

% find closest point one by one
counter      = 0;
flag         = 1;
indx         = zeros (size (X, 1), 2);
tcounter     = 1; % tree counter for tree-by-tree adding
while ~isempty (dplen) && (flag == 1)
    if strfind (options, '-w')
        if mod (counter, 500) == 0,
            waitbar (counter / lenX, HW);
        end
    end
    % find which tree is closest to a point
    distree  = ones (lent, 1);
    for ward = 1 : lent,
        if isempty (dplen {ward}),
            distree(ward) = NaN;
        else
            distree(ward) = min (dplen{ward}, [], 2);
        end
    end
    if sum (isnan (distree)) == lent, % NaN: distance is bigger than tthr
        flag = 0;
    else
        if strfind (options, '-i'), % tree-by-tree addition
            tcounter = tcounter + 1;
            if tcounter > lent
                tcounter = 1;
            end
        else
            % choose the tree that has smallest error:
            [~, tcounter] = min (distree);
        end
        % choose closest point:
        % iopen, index in Open points of vicinity:
        [~, iopen] = min (dplen{tcounter}, [], 2);
        % itree, index in tree:
        itree = ITREE{tcounter}(iopen);
        % update vicinity distance:
        tthr {tcounter} = max ( ...
            tthr{tcounter}, ...
            thr + root_dist{tcounter}(inX{tcounter}(iopen)));
        % update adjacency matrix dA:
        msttrees {tcounter}.dA (end + 1,   itree) = 1;
        msttrees {tcounter}.dA (  itree, end + 1) = 0;
        N {tcounter} = N {tcounter} + 1; % update number of nodes in tree
        % calculate the actual distance of the point to its closest
        % partner in the tree (itree)
        dis = sqrt( ...
            (X (inX{tcounter}(iopen)) - msttrees{tcounter}.X(itree)).^2 + ...
            (Y (inX{tcounter}(iopen)) - msttrees{tcounter}.Y(itree)).^2 + ...
            (Z (inX{tcounter}(iopen)) - msttrees{tcounter}.Z(itree)).^2);
        % don't allow to go beyond the threshold distance:
        dis (dis > thr) = NaN;
        dis (dis > plthr) = NaN;
        % and add this to the path length of that point (itree) to get
        % the total path length to the new point:
        plen_new    = plen {tcounter} (itree) + dis;
        plen{tcounter} = [plen{tcounter}; plen_new];
        if strfind (options, '-p')  % relative path length cost
            disroot = sqrt( ...
                (X (inX {tcounter} (iopen)) - msttrees {tcounter}.X (1)).^2 + ...
                (Y (inX {tcounter} (iopen)) - msttrees {tcounter}.Y (1)).^2 + ...
                (Z (inX {tcounter} (iopen)) - msttrees {tcounter}.Z (1)).^2);            
            plencost {tcounter} = [plencost{tcounter}; ...
                plen_new-disroot];
        else
            plencost {tcounter} = [plencost{tcounter}; ...
                plen_new];
        end        
        % update node coordinates in tree
        msttrees{tcounter}.X = [msttrees{tcounter}.X; ...
            X(inX {tcounter} (iopen))];
        msttrees{tcounter}.Y = [msttrees{tcounter}.Y; ...
            Y(inX {tcounter} (iopen))];
        msttrees{tcounter}.Z = [msttrees{tcounter}.Z; ...
            Z(inX {tcounter} (iopen))];
        msttrees{tcounter}.D = [msttrees{tcounter}.D; ...
            1];
        if ~isempty (L)
            iL = L (inX {tcounter} (iopen));
            msttrees {tcounter}.R = [msttrees{tcounter}.R;  ...
                iL];
        else
            msttrees {tcounter}.R = [msttrees{tcounter}.R; ...
                1];
        end
        % remember which node came from where:
        indx (inX {tcounter} (iopen), :) = [tcounter ...
            length(msttrees {tcounter}.X)];
        if ~isempty (DIST),
            % move node index of DIST matrix from open points to tree:
            iDIST  {tcounter} = [iDIST{tcounter} iDISTP{tcounter}(iopen)];
            iDISTP {tcounter} (iopen) = [];
        end
        % eliminate point in other trees:
        for mation = [(1 : tcounter - 1) (tcounter + 1 : lent)],
            iiopen  = find (inX {mation}    == inX {tcounter} (iopen));
            dplen  {mation} (iiopen)  = [];
            inX    {mation} (iiopen)  = [];
            ITREE  {mation} (iiopen)  = [];
            iiiopen = find (irdist {mation} == inX {tcounter} (iopen));
            irdist {mation} (iiiopen) = [];
            rdist  {mation} (iiiopen) = [];
            if iiiopen <= avic {mation}
                avic {mation} = avic {mation} - 1;
            end;
            if ~isempty (DIST),
                % get rid of indices in DIST of open nodes in all other
                % trees
                iDISTP {mation} (iiopen) = [];
            end
        end
        % get rid of point in open points in vicinity
        dplen{tcounter} (iopen) = [];
        inX   {tcounter} (iopen) = [];
        ITREE{tcounter} (iopen) = [];
        % compare point to dplen to point which is now in the tree
        if ~isempty (dplen{tcounter}), % update in current vicinity
            dis = (sqrt ( ...
                (X (inX {tcounter}) - msttrees {tcounter}.X (end)).^2 + ...
                (Y (inX {tcounter}) - msttrees {tcounter}.Y (end)).^2 + ...
                (Z (inX {tcounter}) - msttrees {tcounter}.Z (end)).^2));
            dis (dis > thr) = NaN;
            dis((dis+plencost{tcounter}(end)) > plthr) = NaN;   % NEW ... if path length of an open point over this point is higher than plthr, ignore it.
            if ~isempty (DIST), % add DISTance matrix factor to Error
                [dplen{tcounter}, idplen] = min ( ...
                    [dplen{tcounter}; ...
                    (dis + bf * plencost {tcounter} (end) + ...
                    Derr * (1 - DIST ( ...
                    iDISTP {tcounter}, ...
                    iDIST {tcounter} (end))))'], ...
                    [], 1);
            else
                [dplen{tcounter}, idplen] = min ( ...
                    [dplen{tcounter}; ...
                    (dis + bf * plencost {tcounter} (end))'], ...
                    [], 1);
            end
            ITREE{tcounter} (idplen == 2) = N {tcounter}; % last added point
            if strfind (options, '-b')
%                 if sum (msttrees {tcounter}.dA (:, itree)) > 1
                    % non-branch points:
                    iCT    = find (sum (msttrees{tcounter}.dA, 1) < 2);
                    if ~isempty (maxT), % there is a max number of stems at root
                        if sum (msttrees{ward}.dA (:, 1), 1) > maxT - 1
                            iCT (iCT == 1) = [];
                        else
                            iCT = unique ([1 iCT]);
                        end
                    end
                    inewbp = find (ITREE{tcounter} == itree);
                    if ~isempty (inewbp),
                        for tetete = 1 : length (inewbp),
                            dis = (sqrt( ...
                                (X (inX  {tcounter} (inewbp (tetete))) - ...
                                msttrees {tcounter}.X (iCT)).^2 + ...
                                (Y (inX  {tcounter} (inewbp (tetete))) - ...
                                msttrees {tcounter}.Y (iCT)).^2 + ...
                                (Z (inX  {tcounter} (inewbp (tetete))) - ...
                                msttrees {tcounter}.Z (iCT)).^2));
                            dis (dis > thr) = NaN;
                            dis (dis+plencost{tcounter}(iCT) > plthr) = NaN;
                            if ~isempty (DIST),
                                [d1 id1] = min ( ...
                                    dis + bf * plencost {tcounter} (iCT) + ...
                                    Derr * (1 - DIST ( ...
                                    iDISTP{tcounter}(inewbp (tetete)), ...
                                    iDIST{tcounter}(iCT) ))', ...
                                    [], 1);
                            else
                                [d1 id1] = min ( ...
                                    dis + bf * plencost {tcounter} (iCT), ...
                                    [], 1);
                            end
                            dplen{tcounter}(inewbp (tetete)) =        d1;
                            ITREE {tcounter} (inewbp (tetete)) = iCT (id1);
                        end
                    end
%                 end
            end
        end

        if ~isempty (L)
            % get rid of all points of one label (axon) if # of synapses full:
            if length (find (msttrees {tcounter}.R == iL)) >= nsyn,
                % first in vicinity and then
                iiL = find (L (inX {tcounter}) == iL);
                % (iiL is index in inX{tcounter} which is part of label to
                % be deleted)
                dplen {tcounter} (iiL) = [];
                inX   {tcounter} (iiL) = [];
                ITREE {tcounter} (iiL) = [];

%                 if ~isempty (avic {tcounter} + 1 : length (rdist{tcounter}))
                    iiL = find (L (irdist {tcounter} (avic{tcounter}+1 : end)) == iL);
                    % (iiL is index in inX{tcounter} which is part of label to
                    % be deleted)
                    rdist  {tcounter} (iiL+avic{tcounter}) = [];
                    irdist {tcounter} (iiL+avic{tcounter}) = [];
%                 end                
%                 if ~isempty (avic {tcounter} + 1 : vic)
%                     iiL = find (L (irdist {tcounter} (avic {tcounter} + ...
%                         1 : vic)) == iL);
%                     % (iiL is index in inX{tcounter} which is part of label to
%                     % be deleted)
%                     rdist  {tcounter} (iiL) = [];
%                     irdist {tcounter} (iiL) = [];
%                 end
            end
        end
        % update vicinity
        vic = sum (rdist {tcounter} < tthr {tcounter});        
        % update dplen etc... according to new vicinity
        if vic > avic {tcounter},
            % new points in vicinity:
            indo = irdist {tcounter} (avic {tcounter} + 1 : vic);
            % number of new points:
            leno = length (indo);
            % repeat the old story with all new points:
            if strfind (options, '-b')
                % non-branch points:
                iCT = find (sum (msttrees{tcounter}.dA, 1) < 2);
                if ~isempty (maxT), % there is a max number of stems at root
                    if sum (msttrees{ward}.dA (:, 1), 1) > maxT - 1
                        iCT (iCT == 1) = [];
                    else
                        iCT = unique ([1 iCT]);
                    end
                end
                dis = sqrt ( ...
                    (repmat (X (indo)', length (iCT), 1) - ...
                    repmat  (msttrees {tcounter}.X (iCT), 1, leno)).^2 + ...
                    (repmat (Y (indo)', length (iCT), 1) - ...
                    repmat  (msttrees {tcounter}.Y (iCT), 1, leno)).^2 + ...
                    (repmat (Z (indo)', length (iCT), 1) - ...
                    repmat  (msttrees {tcounter}.Z (iCT), 1, leno)).^2);
                dis (dis > thr) = NaN;
                dis(dis > plthr) = NaN;
                if ~isempty (DIST),
                    [d1 id1] = min ( ...
                        dis + bf * repmat (plencost {tcounter} (iCT), 1, leno) + ...
                        Derr * (1 - DIST ( ...
                        sub2ind (size (DIST), ...
                        repmat (              indo, length (iCT), 1), ...
                        repmat (iDIST {tcounter} (iCT),         leno, 1)'))), ...
                        [], 1);
                else
                    [d1 id1] = min ( ...
                        dis + bf * repmat (plencost {tcounter} (iCT), 1, leno), ...
                        [], 1);
                end
                id1 = iCT (id1);
            else
                dis = sqrt ( ...
                    (repmat (X (indo)', N {tcounter}, 1) - ...
                    repmat  (msttrees {tcounter}.X, 1, leno)).^2 + ...
                    (repmat (Y (indo)', N {tcounter}, 1) - ...
                    repmat  (msttrees {tcounter}.Y, 1, leno)).^2 + ...
                    (repmat (Z (indo)', N {tcounter}, 1) - ...
                    repmat  (msttrees {tcounter}.Z, 1, leno)).^2);
                dis (dis > thr) = NaN;
                if ~isempty (DIST),
                    [d1 id1] = min ( ...
                        dis + bf * repmat (plencost {tcounter}, 1, leno) + ...
                        Derr * (1 - DIST ( ...
                        sub2ind (size (DIST), ...
                        repmat (        indo, N {tcounter}, 1), ...
                        repmat (iDIST {tcounter},     leno, 1)'))), ...
                        [], 1);
                else
                    [d1 id1] = min ( ...
                        dis + bf * repmat (plencost {tcounter}, 1, leno), ...
                        [], 1);
                end
            end
            dplen {tcounter} = [dplen{tcounter}  d1];
            ITREE {tcounter} = [ITREE{tcounter} id1];
            inX   {tcounter} = [inX{tcounter}  indo];
            if ~isempty (DIST),
                iDISTP {tcounter} = [iDISTP{tcounter} indo+TSUM];
            end
            if strfind (options, '-s'),
                plot3 (X (indo), Y (indo), Z (indo), 'g.');
            end
            avic {tcounter} = vic;
        end
        if strfind (options, '-s'),
            set (HP {tcounter}, 'visible', 'off');
            HP {tcounter} = plot_tree (msttrees {tcounter}, colors (tcounter, :));
            drawnow;
            pause (0.5);
        end
        if strfind (options, '-t'),
            timetrees{tcounter}{end+1} = msttrees {tcounter};
        end
        % indicates that a point has been added in at least one tree
        flag    =           1;
        counter = counter + 1;
    end
end
if strfind (options, '-w'),
    close (HW);
end
if strfind (options, '-t'),
    msttrees = timetrees;
end
if (nargout > 0),
    if lent == 1,
        tree = msttrees {1};
    else
        tree = msttrees;
    end
else
    trees = [trees msttrees];
end

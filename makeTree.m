function [tree,usedpts,availablepts,trueusedpts] = makeTree(availablepts,NodesPerCell,NodesStDev,Somata,bf,pointdens,thicknesses,mode,shapeprofile)
distrtheselayers = 2:4;

% if numel(fac_newPs) ~= distrtheselayers(end)
    fac_newPs = ones(1,distrtheselayers(end));
% end
if nargin < 14
    trueusedpts = {};
end


    NodesN      = round(NodesPerCell + NodesStDev.*randn(1,1));

    % N{1}       = round(normrnd(0.12,0.041) * NodesN); %0.14
    Partrand = 0.6;
    Nrand = cell(5,1);
    %     if isempty(trueusedpts)
    switch mode
        case 1 %AAV rat
            percent = [0 0.20 0.80];
        case 2 %RV rat
            percent = [0.25,0.20,0.55];
        case 3 %AAV mouse
            percent = [0 0.4 0.6];
        case 4 %RV mouse
            percent = [0.25,0.35,0.40];  % not tested yet
    end
    Nrand{2} = percent(1)*NodesN*Partrand;%sum(randist(2:round((numel(randist)-1)/3)))*NodesN*Partrand;
    Nrand{3} = percent(2)*NodesN*Partrand;%sum(randist(round((numel(randist)-1)/3)+1:round((numel(randist)-1)/3*2)))*NodesN*Partrand;
    Nrand{4} = percent(3)*NodesN*Partrand;%sum(randist(round((numel(randist)-1)/3)*2+1:numel(randist)))*NodesN*Partrand;
    %     end
    Ndir{1} = min(3,max(0,round(normrnd(1,0*0.4)))); %0.14  % fixed number in GCL
    if mode <= 2 % -> rat
        Ndir{2} = 0.6*(NodesN*(1-Partrand)-Ndir{1}); %* Nrand{3}/(Nrand{4}+Nrand{3}); % number of directed points in IML dependent on portion of random points in MML compared to OML
        %     Ndir{2} = (NodesN*(1-Partrand)-Ndir{1}) * 0.5 * Nrand{4}/(Nrand{4}+Nrand{3}) +Ndir{2};  % number of directed points in MML dependent on portion of random points in OML compared to MML but points are distributed over IML and MML
        Ndir{3} = 0.15*(NodesN*(1-Partrand)-Ndir{1}); %* 0 * Nrand{4}/(Nrand{4}+Nrand{3}); % number of directed points in MML dependent on portion of random points in OML compared to MML but points are distributed over IML and MML
        Ndir{4} = 0.25*(NodesN*(1-Partrand)-Ndir{1}); %* 0.5 * Nrand{4}/(Nrand{4}+Nrand{3}); % number of directed points in MML dependent on portion of random points in OML compared to MML but points are distributed over IML and MML
    else
        Ndir{2} = 0.3*(NodesN*(1-Partrand)-Ndir{1}); %* Nrand{3}/(Nrand{4}+Nrand{3}); % number of directed points in IML dependent on portion of random points in MML compared to OML
        %     Ndir{2} = (NodesN*(1-Partrand)-Ndir{1}) * 0.5 * Nrand{4}/(Nrand{4}+Nrand{3}) +Ndir{2};  % number of directed points in MML dependent on portion of random points in OML compared to MML but points are distributed over IML and MML
        Ndir{3} = 0.55*(NodesN*(1-Partrand)-Ndir{1}); %* 0 * Nrand{4}/(Nrand{4}+Nrand{3}); % number of directed points in MML dependent on portion of random points in OML compared to MML but points are distributed over IML and MML
        Ndir{4} = 0.15*(NodesN*(1-Partrand)-Ndir{1}); %* 0.5 * Nrand{4}/(Nrand{4}+Nrand{3}); % number of directed points in MML dependent on portion of random points in OML compared to MML but points are distributed over IML and MML
    end
    

pts = cell(4,1);
    % Choose points randomly from each layer
    % pts{1}     = availablepts{1}(randperm(size(availablepts{1}, 1), N{1}), :);
    str = {'','IML','MML','OML'};
    for f = distrtheselayers % no pts in GCL
        ind = randperm(size(availablepts{f}, 1), min(size(availablepts{f}, 1),round(Nrand{f})));
        pts{f} = availablepts{f}(ind,:);
            usedpts{f} = availablepts{f}(ind,:);
  
        availablepts{f}(ind,:) = [];

    end
    pts{4} = cat(1,pts{4},giveOMLpoints(round(Ndir{4}),pointdens{5},Somata,thicknesses,pts{4},shapeprofile));  % these are the directed OML points
    
    pts{3} = cat(1,pts{3},giveMMLpoints(round(Ndir{3}),pointdens{5},Somata,thicknesses,pts{4},shapeprofile));  % these are the directed MML points
    pts{2} = cat(1,pts{2},giveIMLpoints(round(Ndir{2}),pointdens{5},Somata,thicknesses,cat(1,pts{3:4}),shapeprofile));  % these are the directed IML points
    pts{1} = giveGCLpoints3(Ndir{1},pointdens{1},Somata,thicknesses,cat(1,pts{2:4}),shapeprofile);
   

all_pts         = double([pts{3};pts{4}]);
all_startpts         = cat(1,Somata,pts{1}); % ,pts{2}
    maxconn = 500;
    maxpl = 390;

%start with GCL (allow multifurcation there
starttree            = xxMST_tree(1,all_startpts(:,1), ...
    all_startpts(:,2),all_startpts(:,3),...
    [], [], bf, maxconn,[], [], maxpl, '-b');   %nstems  % 0.5 small bf to not have to many branches but to have more curved dend

%add IML
starttree            = xxMST_tree({starttree},pts{2}(:,1), ...
    pts{2}(:,2),pts{2}(:,3),...
    [], [], bf, maxconn,[], [], maxpl, '-b');    %nstems

tree            = xxMST_tree({starttree},all_pts(:,1), ...
    all_pts(:,2),all_pts(:,3),...
    [], [], bf, maxconn,[], [], maxpl, '-b');    %nstems % bf 0.8500000000000000000000000000

display(sprintf('Nodes %d, tree %d',round(sum(cat(1,Nrand{:}))+sum(cat(1,Ndir{:}))),numel(tree.Y)))


for f = distrtheselayers
    availablepts{f} = cat(1,availablepts{f},setdiff(pts{f},[tree.X,tree.Y,tree.Z],'rows'));
    if fac_newPs(f) ~= 1
        usedpts{f} = intersect(usedpts{f},[otree.X,otree.Y,otree.Z],'rows');
    end
end
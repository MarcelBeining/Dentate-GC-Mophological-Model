%% parameters 
targetfolder = pwd;
animal = questdlg('Make mouse or rat DGCs?','What?','Mouse','Rat','Cancel','Mouse');
if isempty(animal) || strcmp(animal,'Cancel')
    return
end
% animal = 'mouse'   % rat
skipRV = 0;
N_cells = 15;

analyzesubsteps = 0;
plotit = 0;
zfac = 1.8;  % 3 factor of growth
xfac = zfac*1.1;
yfac = zfac*1.25;
Suprapyramidal = 1;

cone_alpha = 33;
cone_beta = 22;
bf              = 0.9;

dSGCL = 60/zfac;
if strcmpi(animal,'mouse')
    thresh = 20; %µm or below getting pruned
    STSthresh = 50; %only for analysis
    dGCL = 80/zfac;  % dicke gcl layer
    dML = round(188/zfac /5)*5 ; % dicke ML layer approximated from SH07
    thicknesses = round([dSGCL,dGCL,dML/5,dML*2/5,dML*2/5]);
    NodesPerCell    = 45; %%%%% achtung NOT VALIDATED YET
else
    thresh = 35; %µm or below getting pruned
    STSthresh = 100; %only for analysis
    dGCL = 80/zfac;  % dicke gcl layer
    dML = round(265/zfac / 3)*3; % dicke ML layer
    NodesPerCell    = 25;
    thicknesses = round([dSGCL,dGCL,dML/3,dML/3,dML/3]);
end
NodesStDev      = 0;%3.5;



pointcol=colorme({'light blue','pink','violett','orange'});
DGsize = [350 350 round(dML+dGCL+dSGCL)];
layerZ = cumsum(thicknesses);
%% point densities and creation
CsRV = [];
CsAAV = [];
% this is the GCL distribution
Pointdens{1} = (1:100).^2;
Pointdens{1} = cumsum(Pointdens{1}/sum(Pointdens{1}));
% this is the IML random point distribution
Pointdens{2} = ones(1,101);%((0:100) - 10).^2;   % minus defines where midpoint is
Pointdens{2} = Pointdens{2}/max(Pointdens{2}) + 0.3 ;  % add also some more points in the midth
Pointdens{2} = cumsum(Pointdens{2} / sum(Pointdens{2}));
% this is the MML random & directed point distribution
Pointdens{3} = ((0:100) - 20).^2;   % minus defines where midpoint is
Pointdens{3} = Pointdens{3}/max(Pointdens{3}) + 0.1;  % add also some more points in the midth
Pointdens{3} = cumsum(Pointdens{3} / sum(Pointdens{3}));
% this is the OML random point distribution
Pointdens{4} = ((0:100) - 20).^2;   % minus defines where midpoint is
tmpdens = fliplr((0:100))+50;
tmpdens(tmpdens>100) = 100;
tmpdens = tmpdens/max(tmpdens);
Pointdens{4} = Pointdens{4}/max(Pointdens{4}) + 0.1;  % add also some more points in the midth
Pointdens{4} = Pointdens{4}/max(Pointdens{4});  % 
Pointdens{4} = Pointdens{4} .* tmpdens;
Pointdens{4} = cumsum(Pointdens{4} / sum(Pointdens{4}));
Pointdens{5} = Pointdens{2};
Pointdens{5} = ((0:100) - 100).^2;   % minus defines where midpoint is
Pointdens{5} = Pointdens{5}/max(Pointdens{5});  % add also some more points in the midth
Pointdens{5} = cumsum(Pointdens{5} / sum(Pointdens{5}));

if plotit
    figure,hold all
    title('Point distributions')
    xlabel('layers')
    plot(1:100,Pointdens{1})%/sum(Pointdens{1}))
    plot(101:201,Pointdens{2})%/sum(Pointdens{2}))
    plot(202:302,Pointdens{3})%/sum(Pointdens{3}))
    plot(303:403,Pointdens{4})%/sum(Pointdens{4}))
    plot(303:403,Pointdens{5})%/sum(Pointdens{5}))101:201
end

tree = cell(1,N_cells);
pts = cell(5,1);
M = false(DGsize(1),DGsize(2),DGsize(3));
TotalPoints = 10000;
Nodes{1}   = round(0.25 * TotalPoints);
Nodes{2}   = round(0.14 * TotalPoints);
Nodes{3}   = round(0.12 * TotalPoints);
Nodes{4}   = round(0.49 * TotalPoints);

% figure,hold all
for l = 1:4
    [~, z] = find(diff(cat(2,zeros(Nodes{l},1),repmat(rand(Nodes{l},1),1,numel(Pointdens{l})) <= repmat(Pointdens{l},Nodes{l},1)),1,2));
    z = z/100*thicknesses(l+1)+layerZ(l);  % scale z values to thickness of layer and add thickness of layers below
    pts{l} = cat(2,rand(Nodes{l},1)*DGsize(1),rand(Nodes{l},1)*DGsize(2),z);
    %     plot3(pts{l}(:,1),pts{l}(:,2),pts{l}(:,3),'x')
end
blade = {'infra','supra'};



%% AAV young

conestruct = cell(N_cells,1);
Somata = conestruct;
firsttree = conestruct;
tree = conestruct;
for t=1:N_cells
    if strcmpi(animal,'mouse')
%         Somata{t} = [round(DGsize(1)/2) round(DGsize(2)/2) min(round(normrnd(0.25,0.05)*(layerZ(2)-layerZ(1)))+layerZ(1),layerZ(2)-1)]; % changed to 80µm GCL...
        Somata{t} = [round(DGsize(1)/2) round(DGsize(2)/2) min(round(normrnd(0.625,0.05)*(layerZ(2)-layerZ(1)))+layerZ(1),layerZ(2)-1)];
    else
        Somata{t} = [round(DGsize(1)/2) round(DGsize(2)/2) min(round(normrnd(0.39,0.12)*(layerZ(2)-layerZ(1)))+layerZ(1),layerZ(2)-1)];
    end
    Cone_Height     = DGsize(3)-Somata{t}(3);
    T_radius = sind(cone_alpha/2) * Cone_Height * 2;
    D_radius = sind(cone_beta/2) * Cone_Height * 2;
    % Create scaled elliptical cone
    vlin            = linspace(0,1.99*pi,100);
    ulinorig            = linspace(0,Cone_Height,30);

    [uorig,~]           = meshgrid(ulinorig,vlin);
    ulin = ulinorig;
    addulin = Somata{t}(3)-layerZ(1);  % with this the cone always starts at this cut below somata
    
    indMMLborder = round(find(ulinorig+Somata{t}(3)>= layerZ(4),1,'first'));  %
    indIMLborder = round(find(ulinorig+Somata{t}(3)>= layerZ(3),1,'first'));
    indGCLborder = round(find(ulinorig+Somata{t}(3)>= layerZ(2),1,'first'));
    ulin(round((numel(ulin)-indMMLborder)/2+indMMLborder):end) = ulin(round((numel(ulin)-indMMLborder)/2+indMMLborder));  % cap cone in OML
    
    ulinx = ulin+addulin;
    uliny = ulin+addulin;

    [u,v]           = meshgrid(ulinx,vlin);
    x_cone          = T_radius*sin(v).*(u)/(Cone_Height+addulin);
    [u,v]           = meshgrid(uliny,vlin);
    y_cone          = D_radius*cos(v).*(u)/(Cone_Height+addulin);
    z_cone          = uorig;
    shapeprofile = [T_radius.*(ulinx)/Cone_Height;D_radius.*(uliny)/Cone_Height;ulinorig]';
    
    x_cone = x_cone + Somata{t}(1);
    y_cone = y_cone + Somata{t}(2);
    z_cone = z_cone + Somata{t}(3);
    
    cone            = [reshape(x_cone,[],1),reshape(y_cone,[],1),reshape(z_cone,[],1)];
    nu = numel(ulin);
    nv = numel(vlin);
    conepoly.faces = triangulate_circshape(nu,nv);
    conepoly.vertices = cone;
    gcl.vertices = [0 0 layerZ(1);DGsize(1) 0 layerZ(1);DGsize(1) DGsize(2) layerZ(1);0 DGsize(2) layerZ(1);0 0 layerZ(2);DGsize(1) 0 layerZ(2);DGsize(1) DGsize(2) layerZ(2);0 DGsize(2) layerZ(2)];
    gcl.faces = triangulate_circshape(2,4);
    
    conestruct{t}.shapeprofile = shapeprofile;
    conestruct{t}.gcl = gcl;
    conestruct{t}.conepoly = conepoly;
    
    if plotit && t == 1
        
        figure;hold on
        p=patch(conepoly,'facealpha',0.1,'edgecolor','k','edgealpha',0.2,'facecolor',[1 0.5 0]);
        plot3(cone(:,1),cone(:,2),cone(:,3),'x')
        hold all
        
        hg = patch(gcl,'facealpha',0.2,'edgealpha',0);
        hl(1)=line(gcl.vertices([1 2 3 4 1],1),gcl.vertices([1 2 3 4 1],2),gcl.vertices([1 2 3 4 1],3),'Color','k','LineWidth',2);
        hl(2)=line(gcl.vertices([5 6 7 8 5],1),gcl.vertices([5 6 7 8 5],2),gcl.vertices([5 6 7 8 5],3),'Color','k','LineWidth',2);
        axis off
    end
    %
    for l = 2:4
        yes = inpolyhedron(conepoly,pts{l},'FLIPNORMALS',1);
        availablepts{l} = pts{l}(yes,:);
    end
    
    [firsttree{t},~,~,trueusedpts{t}] = makeTree(availablepts,NodesPerCell,NodesStDev,Somata{t},bf,Pointdens,thicknesses,1+2*strcmpi(animal,'mouse'),conestruct{t}.shapeprofile);
    firsttree{t}.rnames = {'SGCL','GCL','IML','MML','OML','outside'};
    numel(firsttree{t}.Y)
    if plotit && t == 1
        for l = 1:4
            hp(l) = plot3(trueusedpts{t}{l}(:,1),trueusedpts{t}{l}(:,2),trueusedpts{t}{l}(:,3),'.','markersize',16,'Color',pointcol{l});
        end
    end
    CsAAV(t) = sum(C_tree(firsttree{t}));
    tree{t} = preview_tree(firsttree{t},'-j');  % CAUTION! ONLY FOR PREVIEW! NO MODELING
    
    if plotit && t == 1
        view(90,0)
        p=plot_tree(tree{t});
    end
    
    if plotit && t == 1
        %     plot_tree(tree{t})
        tprint(fullfile(targetfolder,'MST_2'),'-HR-png')
        delete(p)
        tprint(fullfile(targetfolder,'MST_1'),'-HR-png')
        p=plot_tree(tree{t});
        tprint(fullfile(targetfolder,'MST_3'),'-HR-png')
        delete(p)
        delete(hp)
        delete(hg)
        delete(hl)
    end
    indsoma = find(strcmp(tree{t}.rnames,'soma'));
    tree{t}.R(tree{t}.Z < layerZ(1) & tree{t}.R~=indsoma) = 1;
    for d = 1:numel(layerZ)-1
        tree{t}.R(tree{t}.Z >= layerZ(d) & tree{t}.Z <= layerZ(d+1) & tree{t}.R~=indsoma) = d+1;
    end
    tree{t}.R(tree{t}.Z > layerZ(5) & tree{t}.R~=indsoma) = 6;
    sum(B_tree(tree{t}) & tree{t}.R == 4)
    
    tree{t}.animal = 'Artificial';
    tree{t}.pyramidal_blade = blade{Suprapyramidal+1};
    tree{t}.lateral_side = 'contra';
    tree{t}.arc = 'Unknown';
    
    % create contours describing the layering
    tree{t}.contours{1,1}.Vertices = interp_border([0,layerZ(1);DGsize(1),layerZ(1)],1,2);
    for l = 2:5
        tree{t}.contours{1,l}.Vertices = interp_border([0,layerZ(l-1);DGsize(1),layerZ(l-1)],1,2);
        tree{t}.contours{1,l}.Vertices = cat(1,tree{t}.contours{1,l}.Vertices,interp_border([DGsize(1),layerZ(l);0,layerZ(l)],1,2));
    end
    % turn tree in Y direction...necessary for analysis script
    tmp = tree{t}.Y;
    tree{t}.Y = tree{t}.Z;
    tree{t}.Z = tmp;

end
%
if analyzesubsteps
    dstruct = Analyze_Morph_regionally(tree,{'Total','IML','MML','OML','outside'},[],'-o-g-c-NPS',[],STSthresh);
    save(sprintf('D:/%s_AAVart_young.mat',animal),'tree','dstruct','Nodes','Pointdens','DGsize','thicknesses','layerZ')
else
    delete(sprintf('D:/%s_AAVart_young.mat',animal))
    save(sprintf('D:/%s_AAVart_young.mat',animal),'tree','Nodes','Pointdens','DGsize','thicknesses','layerZ')
end

%% prune AAV tree
load(sprintf('D:/%s_AAVart_young.mat',animal),'tree')
for t = 1:numel(tree)
    PL = Pvec_tree(tree{t});
    [tree{t},AAVcount(t)] = prune_tree(tree{t},thresh);%,setdiff(1:5,find(strcmp(tree{t}.rnames,'OML'))));%,'-s');
end
if analyzesubsteps
    dstruct = Analyze_Morph_regionally(tree,{'Total','IML','MML','OML','outside'},[],'-o-g-c-NPS',[],STSthresh);
    save(sprintf('D:/%s_AAVart_young_pruned.mat',animal),'tree','dstruct','DGsize','thicknesses','conestruct','Nodes','Pointdens','DGsize','thicknesses','layerZ')
else
    delete(sprintf('D:/%s_AAVart_young_pruned.mat',animal))
    save(sprintf('D:/%s_AAVart_young_pruned.mat',animal),'tree','Nodes','Pointdens','DGsize','thicknesses','layerZ')
end
%% DG growth
% r1 = (layerZ(end)-layerZ(1))* 10 *10000;  % approximated radius of transverse curvature of dg (thickness of dg is about one tenth (betw. 1/10 and 1/5) of curvature radius to SGCL), in crest much much higher!
% r2 = (layerZ(end)-layerZ(1))* 20 *10000;  %approximated radius for longitudianal curvature of dg



DGsize(3) = round(DGsize(3) * zfac);
thicknesses = round(thicknesses * zfac);
layerZ = round(layerZ * zfac);
gcl.vertices = [0 0 layerZ(1);DGsize(1) 0 layerZ(1);DGsize(1) DGsize(2) layerZ(1);0 DGsize(2) layerZ(1);0 0 layerZ(2);DGsize(1) 0 layerZ(2);DGsize(1) DGsize(2) layerZ(2);0 DGsize(2) layerZ(2)];
gcl.faces = triangulate_circshape(2,4);
gcl.vertices = gcl.vertices(:,[1 3 2]);
for t=1:N_cells
    
    conestruct{t}.shapeprofile(:,3) = conestruct{t}.shapeprofile(:,3) * zfac;
    conestruct{t}.gcl.vertices(:,3) = conestruct{t}.gcl.vertices(:,3) * zfac;
    
    AAVconestruct{t} = conestruct{t};
    AAVconestruct{t}.conepoly.vertices(:,1) = (AAVconestruct{t}.conepoly.vertices(:,1)-Somata{t}(1))  .* xfac  + Somata{t}(1); %.* (r1+ AAVconestruct{t}.conepoly.vertices(:,3)*zfac)./(r1+AAVconestruct{t}.conepoly.vertices(:,3)))
    AAVconestruct{t}.conepoly.vertices(:,2) = (AAVconestruct{t}.conepoly.vertices(:,2)-Somata{t}(2))  .* yfac  + Somata{t}(2); %.* (r2+ AAVconestruct{t}.conepoly.vertices(:,3)*zfac)./(r2+AAVconestruct{t}.conepoly.vertices(:,3)))
    AAVconestruct{t}.conepoly.vertices(:,3) = AAVconestruct{t}.conepoly.vertices(:,3) * zfac;
    
    conestruct{t}.conepoly.vertices(:,3) = conestruct{t}.conepoly.vertices(:,3) * zfac;
    
    tree{t} = tran_tree(tree{t},-[Somata{t}(1),0,Somata{t}(2)]); %layerZ(2)
    % keep in mind that Y is the new Z in the old trees
    tree{t}.X = tree{t}.X  .* xfac;  % .* (r1+tree{t}.Y*zfac)./(r1+tree{t}.Y)    this is the same spread calculation as the spread above + an extra xyfac elongation
    tree{t}.Z = tree{t}.Z  .* yfac;   %.* (r2+tree{t}.Y*zfac)./(r2+tree{t}.Y)      " + an extra xyfac elongation
    tree{t}.Y = tree{t}.Y * zfac;
    tree{t} = tran_tree(tree{t},[Somata{t}(1),0,Somata{t}(2)]); %newlayerZ(2)
    %%%%%%%%%%%%%
%     figure(666),hold on
%     plot_tree(tree{t})
%     tree{t} = restrain_tree(tree{t},maxpl);
%     plot_tree(tree{t},[1 0 0])
%     pause
%     clf
    %%%%%%%%%%%%%
    % create contours describing the layering
    tree{t}.contours{1,1}.Vertices = interp_border([0,layerZ(1);DGsize(1),layerZ(1)],1,2);
    for l = 2:5
        tree{t}.contours{1,l}.Vertices = interp_border([0,layerZ(l-1);DGsize(1),layerZ(l-1)],1,2);
        tree{t}.contours{1,l}.Vertices = cat(1,tree{t}.contours{1,l}.Vertices,interp_border([DGsize(1),layerZ(l);0,layerZ(l)],1,2));
    end
    if strcmpi(animal,'mouse')
        tree{t} = preview_tree(tree{t},'-m');  % achtung, jitter darf nach region zuordnung kommen!
    else
        tree{t} = preview_tree(tree{t},'-r');  % achtung, jitter darf nach region zuordnung kommen!
    end
end

dstruct = Analyze_Morph_regionally(tree,{'Total','IML','MML','OML','outside'},[],'-o-g-c-NPS',[],STSthresh);
save(sprintf('D:/%s_AAVart_old_pruned.mat',animal),'tree','dstruct','firsttree','AAVconestruct','Nodes','Pointdens','DGsize','thicknesses','layerZ')
save_tree(tree,fullfile(targetfolder,sprintf('%s_AAVart_old_pruned.mtr',animal)))

%% RV tree
if skipRV
    return
end

load(sprintf('D:/%s_AAVart_old_pruned.mat',animal))
for l = 1:4
    [~, z] = find(diff(cat(2,zeros(Nodes{l},1),repmat(rand(Nodes{l},1),1,numel(Pointdens{l})) <= repmat(Pointdens{l},Nodes{l},1)),1,2));
    z = z/100*thicknesses(l+1)+layerZ(l);  % scale z values to thickness of layer and add thickness of layers below
    pts{l} = cat(2,rand(Nodes{l},1)*DGsize(1),rand(Nodes{l},1)*DGsize(2),z);
    %     plot3(pts{l}(:,1),pts{l}(:,2),pts{l}(:,3),'x')
end
for t = 1:N_cells
    Somata{t} = [round(DGsize(1)/2) round(DGsize(2)/2) min(round(normrnd(0.17,0.12)*(layerZ(2)-layerZ(1)))+layerZ(1),layerZ(2)-1)];
    
    Cone_Height     = DGsize(3)-Somata{t}(3);
    T_radius = sind(cone_alpha/2) * Cone_Height * 2;
    fprintf('T_radius %g\n',T_radius)
    D_radius = sind(cone_beta/2) * Cone_Height * 2;
    % Create scaled elliptical cone
    vlin            = linspace(0,1.99*pi,100);
    ulinorig            = linspace(0,Cone_Height,30);

    [uorig,~]           = meshgrid(ulinorig,vlin);
    ulin = ulinorig;
    addulin = Somata{t}(3)-layerZ(1);  % with this the cone always starts at this cut below somata
    ulinx = ulin+addulin;
    uliny = ulin+addulin;

    [u,v]           = meshgrid(ulinx,vlin);
    x_cone          = T_radius*sin(v).*(u)/(Cone_Height+addulin);
    [u,v]           = meshgrid(uliny,vlin);
    y_cone          = D_radius*cos(v).*(u)/(Cone_Height+addulin);
    z_cone          = uorig;
    shapeprofile = [T_radius.*(ulinx)/Cone_Height;D_radius.*(uliny)/Cone_Height;ulinorig]';
    
    x_cone = x_cone + Somata{t}(1);
    y_cone = y_cone + Somata{t}(2);
    z_cone = z_cone + Somata{t}(3);
    
    cone            = [reshape(x_cone,[],1),reshape(y_cone,[],1),reshape(z_cone,[],1)];
    nu = numel(ulin);
    nv = numel(vlin);
    conepoly.faces = triangulate_circshape(nu,nv);
    conepoly.vertices = cone;
    gcl.vertices = [0 0 layerZ(1);DGsize(1) 0 layerZ(1);DGsize(1) DGsize(2) layerZ(1);0 DGsize(2) layerZ(1);0 0 layerZ(2);DGsize(1) 0 layerZ(2);DGsize(1) DGsize(2) layerZ(2);0 DGsize(2) layerZ(2)];
    gcl.faces = triangulate_circshape(2,4);
    
    conestruct{t}.shapeprofile = shapeprofile;
    conestruct{t}.gcl = gcl;
    conestruct{t}.conepoly = conepoly;
    
    if plotit && t == 1
        
        figure;hold on
        hold all
        
        hg = patch(gcl,'facealpha',0.2,'edgealpha',0);
        hl(1)=line(gcl.vertices([1 2 3 4 1],1),gcl.vertices([1 2 3 4 1],2),gcl.vertices([1 2 3 4 1],3),'Color','k','LineWidth',2);
        hl(2)=line(gcl.vertices([5 6 7 8 5],1),gcl.vertices([5 6 7 8 5],2),gcl.vertices([5 6 7 8 5],3),'Color','k','LineWidth',2);
        axis off
    end
    %
    for l = 2:4
        yes = inpolyhedron(conepoly,pts{l},'FLIPNORMALS',1);
        %         sum(yes)
%         h=plot3(pts{l}(yes,1),pts{l}(yes,2),pts{l}(yes,3),'x');
        availablepts{l} = pts{l}(yes,:);
    end
    
    [tree{t}] = makeTree(availablepts,NodesPerCell,NodesStDev,Somata{t},bf,Pointdens,thicknesses,2+2*strcmpi(animal,'mouse'),conestruct{t}.shapeprofile);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tree{t} = restrain_tree(tree{t},maxpl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tree{t}.rnames = {'SGCL','GCL','IML','MML','OML','outside'};    
    CsRV(t) = sum(C_tree(tree{t}));
    if strcmpi(animal,'mouse')
        tree{t} = preview_tree(tree{t},'-m-j');  % achtung, jitter darf nach region zuordnung kommen!
    else
        tree{t} = preview_tree(tree{t},'-r-j');  % achtung, jitter darf nach region zuordnung kommen!
    end
    if plotit && t == 1
        patch(conestruct{t}.gcl,'facealpha',0.2,'edgealpha',0)
        hl(1)=line(conestruct{t}.gcl.vertices([1 2 3 4 1],1),conestruct{t}.gcl.vertices([1 2 3 4 1],2),conestruct{t}.gcl.vertices([1 2 3 4 1],3),'Color','k','LineWidth',2);
        hl(2)=line(conestruct{t}.gcl.vertices([5 6 7 8 5],1),conestruct{t}.gcl.vertices([5 6 7 8 5],2),conestruct{t}.gcl.vertices([5 6 7 8 5],3),'Color','k','LineWidth',2);
        plot_tree(tree{t})
        tprint(fullfile(targetfolder,'MST_4'),'-HR-png')
    end
    
    indsoma = find(strcmp(tree{t}.rnames,'soma'));
    tree{t}.R(tree{t}.Z < layerZ(1) & tree{t}.R~=indsoma) = 1;
    for d = 1:numel(layerZ)-1
        tree{t}.R(tree{t}.Z >= layerZ(d) & tree{t}.Z <= layerZ(d+1) & tree{t}.R~=indsoma) = d+1;
    end
    tree{t}.R(tree{t}.Z > layerZ(5) & tree{t}.R~=indsoma) = 6;
    
    
    tree{t}.animal = 'Artificial';
    tree{t}.pyramidal_blade = blade{Suprapyramidal+1};
    tree{t}.lateral_side = 'contra';
    tree{t}.arc = 'Unknown';
    
    % create contours describing the layering
    tree{t}.contours{1,1}.Vertices = interp_border([0,layerZ(1);DGsize(1),layerZ(1)],1,2);
    for l = 2:5
        tree{t}.contours{1,l}.Vertices = interp_border([0,layerZ(l-1);DGsize(1),layerZ(l-1)],1,2);
        tree{t}.contours{1,l}.Vertices = cat(1,tree{t}.contours{1,l}.Vertices,interp_border([DGsize(1),layerZ(l);0,layerZ(l)],1,2));
    end
    % turn tree in Y direction...necessary for analysis script
    tmp = tree{t}.Y;
    tree{t}.Y = tree{t}.Z;
    tree{t}.Z = tmp;
end
if analyzesubsteps
    dstruct = Analyze_Morph_regionally(tree,{'Total','IML','MML','OML','outside'},[],'-o-g-c-NPS',[],STSthresh);
    save(sprintf('D:/%s_RVart.mat',animal),'tree','dstruct','conestruct')
else
    delete(sprintf('D:/%s_RVart.mat',animal))
    save(sprintf('D:/%s_RVart.mat',animal),'tree','conestruct')
end
%% prune RV  tree
load(sprintf('D:/%s_RVart.mat',animal))

for t = 1:numel(tree)
    PL = Pvec_tree(tree{t});
    [tree{t},RVcount(t)] = prune_tree(tree{t},thresh,[]);%,setdiff(1:5,find(strcmp(tree{t}.rnames,'OML')))); %prune everywhere except OML
end
dstruct = Analyze_Morph_regionally(tree,{'Total','IML','MML','OML','outside'},[],'-o-g-c-NPS',[],STSthresh);
save(sprintf('D:/%s_RVart_pruned.mat',animal),'tree','dstruct','conestruct')
save_tree(tree,fullfile(targetfolder,sprintf('%s_RVart_pruned.mtr',animal)))

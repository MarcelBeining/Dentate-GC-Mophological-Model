function tree = preview_tree(tree,options)
%-s show
%-r rat
%-m mouse

if nargin < 2
    options = '-s';
end

if ~isempty(strfind(options,'-j'))
    
    
    % Define spatial jitter amplitude
    % stde_mean   = 0.33;%.3;%0.25
    % stde_stdev  = 0;%.125;%0.125
    stde_mean2   = 0.4;%.2;%.5;
    stde_stdev2  = 0;%0.1;    %0.05; %0.1
    stde_mean3   = 0.23;%.1;%0.05;%0.15%0.25
    stde_stdev3  = 0;%0.06;    %0.06;%0.125
    
    
    % Jitter tree
    % stde  = normrnd(stde_mean,stde_stdev);
    stde2 = normrnd(stde_mean2,stde_stdev2);
    stde3 = normrnd(stde_mean3,stde_stdev3);
    % figure
    tree = smooth_tree(tree,[],0.5,1,' '); % to reduce strong turnings
    tree = resample_tree(tree,4,'-r-b');  % caution! resampling has huge effect on jittering!!!
    % tree = jitter_tree(tree,stde,100,' ');%,find(tree.R == find(strcmp(tree.rnames,'SGCL')) | tree.R == find(strcmp(tree.rnames,'GCL'))));
    tree = jitter_tree(tree,stde2,10,' ');
    tree = jitter_tree(tree,stde3,4,' ');
end
tree = resample_tree(tree,1,'-r-b');

if ~any(strcmp(tree.rnames,'soma'))
    tree.rnames = [tree.rnames,'soma'];
    R = numel(tree.rnames);
else
    R = find(strcmp(tree.rnames,'soma'));
end

if ~isempty(strfind(options,'-r'))
    tree = quaddiameter_tree(tree,normrnd(0.08,0.008),normrnd(0.5,0.05));%,0.091,0.381);  % der neuere wert ist korrektur für  schwache signalstärke bei 2P für distale dendriten... % die werte sind die gemittelten fits aus den AAV contra (supra+infra) Rekonstruktionen
    tree = soma_tree(tree,normrnd(21.5,1),25,'-b');%,'-r'); % also add region soma, 22 means 11 µm diameter...
    tree.R(Pvec_tree(tree) < 8) = R;

elseif ~isempty(strfind(options,'-m'))
    %     tree = jitter_tree(tree,stde3,3,' ');
    %     tree = resample_tree(tree,1,'-r-b');
    tree = quaddiameter_tree(tree,normrnd(0.100,0.01),normrnd(0.396,0.0396));%0.11,0.65);  % die neuen werte sind die gemittelten fits aus allen SH07 zellen.
    tree = soma_tree(tree,normrnd(20,1),23,'-b');%,'-r'); % also add region soma, 22 means 11 µm diameter...
    tree.R(Pvec_tree(tree) < 4.2) = R; %21/5

else
    tree.D(:) = 2;
end



if ~isempty(strfind(options,'-s'))
    plot_tree(tree)
end

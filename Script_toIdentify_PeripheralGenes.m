clc
clear
close all

dir =       '...';
CCflag =    1;

diseases = 1:70;

% % Setting the variables for the elaborations
% sdir = fullfile('M:\Network Medicine\DA&R\Results_DiaBLe',sprintf('HIPPIE thr%g no self-loops',thrHscore));
% if ~exist(sdir,'dir')
%     mkdir(sdir);
% end

%% Importing data
%Import the interactome 
Interactome = readtable('...\PPI201806_large.txt');
genes(:,1) = Interactome.Var1;
genes(:,2) = Interactome.Var3;
%save the gene names in a variable
gene_names(:,1) = Interactome.Var2;
gene_names(:,2) = Interactome.Var4;

[all_genes,ia,ic] = unique(genes,'stable');
[all_genes_names,ian,icn] = unique(gene_names,'stable');
N = length(all_genes);% total number of nodes in the entire network
clear Interactome

%% Creating the adjacency matrix
adj=zeros(N);
for i=1:N
    ind = find(genes(:,1)==all_genes(i,1));
    ind2 = find(genes(:,2)==all_genes(i,1));
    if ~isempty(ind)
        for j=1:length(ind)
            ind_con = find(all_genes(:,1)==genes(ind(j,1),2));
            adj(i,ind_con) = 1;
            adj(ind_con,i) = 1;
        end
    end
    if ~isempty(ind2)
        for j=1:length(ind2)
            ind_con2 = find(all_genes(:,1)==genes(ind2(j,1),1));
            adj(i,ind_con2) = 1;
            adj(ind_con2,i) = 1;
        end
    end
    clear ind_con ind_con2 ind ind2
end
adj = adj-triu(tril(adj));
%check for isolated nodes and elimination of them
K = sum(adj,2);
adj2 = adj;
adj2(K==0,:) = []; 
adj2(:,K==0) = [];
wholeG = graph(adj2);
all_genes(K==0) = [];
% all_genes_names(K==0) = [];
clear adj adj2 K

%% to identify peripheral nodes considering the whole graph or the giant connected compontent
if CCflag==1
    [bins,binsizes] = conncomp(wholeG);
    idx = binsizes(bins) == max(binsizes);
    G = subgraph(wholeG, idx);
    all_genes = all_genes(idx);
%     all_genes_names = all_genes_names(idx);
else
    G = wholeG;
end

N = length(all_genes);% total number of nodes in the entire network (isolated nodes are excluded)

%% importing the seeds for 70 diseases from a tsv file
disease = tdfread('seeds.tsv');
genes_seed = cell(size(disease.Genes,1),1);
for i=1:70
    ge = regexp(disease.Genes(i,:),'\\(\d+)','tokens');
    for j=1:length(ge)
        genes_seed{i,1}(1,j) = str2num(ge{1,j}{1,1});
    end
end
clear ge
SeedsPeriphGenes.diseases = disease.Diseases;
SeedsPeriphGenes.seeds = genes_seed;

%%
for dd=diseases
    %Chose the disease that you prefere to look for
    seeds0 = genes_seed{dd,1}; 
    [~, locS] = ismember(seeds0,all_genes);
    locS(find(locS==0)) = [];
    seeds = locS;
    PeriphNodes = setdiff(1:numel(all_genes),seeds);
    while numel(seeds0)/numel(PeriphNodes)<10^-1
        numel(seeds0)/numel(PeriphNodes)
        
        AllNeighb = [];
        for ss=1:length(seeds)
            Neighb = neighbors(G,seeds(ss));
            Neighb = setdiff(Neighb,seeds);
            AllNeighb = [AllNeighb;Neighb];
            clear Neighb
        end
        AllNeighbors = unique(AllNeighb,'stable');
        seeds = union(seeds,AllNeighbors);
        PeriphNodes = setdiff(PeriphNodes,AllNeighbors);
        clear AllNeighb AllNeighbors
    end
   periph_genes{dd} = all_genes(PeriphNodes);
   
   clear locS PeriphNodes
     
end

SeedsPeriphGenes.periph_genes = periph_genes';

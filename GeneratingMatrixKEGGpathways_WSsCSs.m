% Script to generate the N x 2M matrix with:
% -> N = genes number 
% -> M = number of significant KEGG pathways for seed genes
% The first half of the GO annotations is related to warm seeds, 
% the other half to cold seeds.

clc
clear

%%
dir =       '...';
savedir =   fullfile(dir,'KEGG Matrices WSs-CSs');
if ~exist(savedir,'dir')
    mkdir(savedir);
end

%% loading data
load('SelectedDiseases_WSsCSs_GOKegg.mat')
load('...\KEGG_pathways.mat')
PWnames = fieldnames(PWgenes);

%% Importing files: Interactome
%Import the interactome 
Interactome = readtable('...\PPI201806_large.txt');
genes(:,1) = Interactome.Var1;
genes(:,2) = Interactome.Var3;
%save the gene names in a variable
gene_names(:,1) = Interactome.Var2;
gene_names(:,2) = Interactome.Var4;

[all_genes,ia,ic] = unique(genes,'stable');
gnames = loadname('...\PPI201806_GenesSymbols.mat');
clear Interactome

TN = length(gnames);

%% Creating the adjacency matrix
adj=zeros(TN);
for i=1:TN
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
adj2(find(K==0),:) = []; adj2(:,find(K==0)) = [];
all_genes(find(K==0)) = [];
gnames(find(K==0)) = [];
clear adj adj2 K NT
M = length(gnames);

%% Importing files: seed sets
%importing the seeds for 70 diseases from a tsv file
disease = tdfread('seeds.tsv');
genes_seed = cell(size(disease.Genes,1),1);
for i=1:70
    ge = regexp(disease.Genes(i,:),'\\(\d+)','tokens');
    for j=1:length(ge)
        genes_seed{i,1}(1,j) = str2num(ge{1,j}{1,1});
    end
end
clear ge

%% Generating the Nx2M matrix
for dd=1:length(SelDiseases.ID)
    
    if ~isempty(SelDiseases.KEGGIndex{dd}) && ~isempty(SelDiseases.PeriphGKEGGIndex{dd})
        project = SelDiseases.ID{dd};
        pathways = [SelDiseases.KEGGIndex{dd}' SelDiseases.PeriphGKEGGIndex{dd}'];
        Matrix = zeros(M,length(pathways));
        for pp=1:length(pathways)
            Matrix(:,pp) = ismember(gnames,PWgenes.(PWnames{pathways(pp)}));
            KEGGpathways{pp} = PWnames{pathways(pp)};
        end
        disgene_seed = genes_seed{eval(project),1};
        seedgenes = intersect(disgene_seed,all_genes);
        clear disgene_seed
        Note = 'The first half of the KEGG pathways is related to WSs, the other half to the CSs';
        
        % saving matlab files
        GenePathMatrix.Matrix = Matrix;
        GenePathMatrix.geneIDs = all_genes;
        GenePathMatrix.geneSymbols = gnames;
        GenePathMatrix.diseasegenes = seedgenes;
        GenePathMatrix.KEGGpathways = KEGGpathways;
        GenePathMatrix.Note = Note;
        save(sprintf('%s\\%s_KEGGpathways_Matrix.mat',savedir,project),'GenePathMatrix')
        
        clear Matrix pathways project KEGGpathways GenePathMatrix seedgenes
    end
    
end

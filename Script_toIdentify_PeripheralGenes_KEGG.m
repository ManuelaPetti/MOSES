
clc
clear

%%
dir =       '...';

%% loading data
load('SelectedDiseases.mat')
load('KEGG_pathways.mat')
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

% [all_genes_names,ian,icn] = unique(gene_names,'stable');
clear Interactome


%% Generating the Nx2M matrix
for dd=1:length(SelDiseases.ID)
    SeedsG = SelDiseases.seed{dd};
    ind_s = ismember(all_genes,SeedsG);
    PeriphG = SelDiseases.periph_genes{dd};
    ind_pg = ismember(all_genes,PeriphG);
    gnamesD = cat(1,gnames(ind_s),gnames(ind_pg));
    
    M = numel(gnamesD);
    if numel(SelDiseases.KEGGIndex{dd})>1
        project = SelDiseases.ID{dd};
        pathways = SelDiseases.KEGGIndex{dd}';
        Matrix = zeros(M,length(pathways));
        for pp=1:length(pathways)
            Matrix(:,pp) = ismember(gnamesD,PWgenes.(PWnames{pathways(pp)}));
        end
        check = sum(Matrix(numel(SeedsG)+1:end,:),2);
        indout = find(check~=0);
        PeriphG(indout) = [];
        NSPG(dd,:) = [numel(SeedsG) numel(PeriphG) numel(pathways) numel(indout)];
        SelDiseases.periph_genes{dd,1} = PeriphG;
    end
    clear SeedsG PeriphG indout check annotations Matrix gnamesD
end

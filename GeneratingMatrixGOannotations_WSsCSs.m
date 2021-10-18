% Script to generate the N x 2M matrix with:
% -> N = genes number 
% -> M = number of significant GO annotations (biological process) for seed genes
% The first half of the GO annotations is related to warm seeds, 
% the other half to cold seeds.

clc
clear

%%
dir =       '...';
savedir =   fullfile(dir,'Gene Ontology Annotations Matrices WSs-CSs');
if ~exist(savedir,'dir')
    mkdir(savedir);
end

propagation =   'no propagation';   %'no propagation' 'getancestors' 'getrelatives'
proj1 =         'BP';               %'BP'-->Aspect','P' in line 77
                                    %'CC'-->Aspect','C' in line 77
                                    %'MF'-->Aspect','F' in line 77
if strcmp(proj1,'BP')
    aspect = 'P';
elseif strcmp(proj1,'CC')
    aspect = 'C';
elseif strcmp(proj1,'MF')
    aspect = 'F';
end


%% loading data
load('SelectedDiseases_WSsCSs_GOKegg.mat')

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

%% Update the Gene Ontology annotations
GO = geneont('live',true);
% Creating the Human Annotation Map
%Open the annotation file found in
%(http://www.geneontology.org/page/download-annotations) selecting only the
%needed characteristics of the annotations, only for molecular functions
%(Aspect=F)
Humann = goannotread('goa_human_valid.gaf','Aspect',aspect,...
    'Fields',{'DB_Object_Symbol','GOid','Evidence','Qualifier'});


ind = [];
%eliminating all the IPI to avoid circularity and the non empty Qualifier
for i=1:size(Humann,1)
    if strcmp(Humann(i).Evidence,'IPI')
        ind = [ind;i];
    end
end
Humann(ind) = []
ind=[];
for i=1:length(Humann)
    if any(Humann(i).Qualifier)
        ind = [ind,i];
    end
end
Humann(ind) = []

%%
%creating a new map that connects to the gene_name all its GO IDs, the keys
%are the gene names, to which can be associated many GO IDs
GOannotations = struct;
for i=1:numel(Humann)
    GOid = Humann(i).GOid;
    %check if this key is part of the map
    if isfield(GOannotations,strcat('GO',num2str(GOid)))
        %Adding the key and the GOid
        GOannotations.(strcat('GO',num2str(GOid))){length(GOannotations.(strcat('GO',num2str(GOid))))+1} = ...
            Humann(i).DB_Object_Symbol;
    else
        GOannotations.(strcat('GO',num2str(GOid))){1} = Humann(i).DB_Object_Symbol;
    end
    clear GOid
end

clear Humann
%%
fldnames = fieldnames(GOannotations);
for ik=1:length(fldnames)
    GOannotations.(fldnames{ik}) = unique(GOannotations.(fldnames{ik}),'stable');
end

%% Generating the Nx2M matrix
for dd=1:length(SelDiseases.ID)
    
    if ~isempty(SelDiseases.GOIndex{dd}) && ~isempty(SelDiseases.PeriphGGOIndex{dd})
        project = SelDiseases.ID{dd};
        seedgenes = SelDiseases.seed{dd};
        annotations = [SelDiseases.GOIndex{dd}' SelDiseases.PeriphGGOIndex{dd}'];
        Matrix = zeros(M,length(annotations));
        for pp=1:length(annotations)
            Matrix(:,pp) = ismember(gnames,GOannotations.(strcat('GO',num2str(annotations(pp)))));
        end
        indseeds = ismember(all_genes,seedgenes);
%         Matrix(indseeds,1:numel(SelDiseases.GOIndex{dd})) = 1;
        seedgenes = intersect(seedgenes,all_genes);
        clear disgene_seed
        Note = 'The first half of the GO annotations is related to WSs, the other half to the CSs';
        
        % saving matlab files
        GeneGOannMatrix.Matrix = Matrix;
        GeneGOannMatrix.geneIDs = all_genes;
        GeneGOannMatrix.geneSymbols = gnames;
        GeneGOannMatrix.diseasegenes = seedgenes;
        GeneGOannMatrix.GOids = annotations;
        GeneGOannMatrix.Note = Note;
        save(sprintf('%s\\%s_GOannotations_Matrix.mat',savedir,project),'GeneGOannMatrix')
     
        clear Matrix annotations project GeneGOannMatrix seedgenes
    end
    
end



clc
clear

%%
dir =       '...';

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
load('SelectedDiseases.mat')

%% Importing files: Interactome
%Import the interactome 
Interactome = readtable('M:\Network Medicine\Data\PPI201806_large\PPI201806_large.txt');
genes(:,1) = Interactome.Var1;
genes(:,2) = Interactome.Var3;
%save the gene names in a variable
gene_names(:,1) = Interactome.Var2;
gene_names(:,2) = Interactome.Var4;

[all_genes,ia,ic] = unique(genes,'stable');
gnames = loadname('M:\Network Medicine\PPI201806_GenesSymbols.mat');

% [all_genes_names,ian,icn] = unique(gene_names,'stable');
clear Interactome

%% Importing files: seed sets & peripheral genes
SeedsPeriphGenes = loadname('...\SeedsPeriphGenes.mat');

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
    SeedsG = SeedsPeriphGenes.seeds{str2num(SelDiseases.ID{dd})};
    ind_s = ismember(all_genes,SeedsG);
    PeriphG = SeedsPeriphGenes.periph_genes{str2num(SelDiseases.ID{dd})};
    ind_pg = ismember(all_genes,PeriphG);
    gnamesD = cat(1,gnames(ind_s),gnames(ind_pg));
    
    M = numel(gnamesD);
    if numel(SelDiseases.GOIndex{dd})>1
        project = SelDiseases.ID{dd};
        annotations = SelDiseases.GOIndex{dd}';
        Matrix = zeros(M,length(annotations));
        for pp=1:length(annotations)
            Matrix(:,pp) = ismember(gnamesD,GOannotations.(strcat('GO',num2str(annotations(pp)))));
        end
        check = sum(Matrix(numel(SeedsG)+1:end,:),2);
        indout = find(check~=0);
        PeriphG(indout) = [];
        NSPG(dd,:) = [numel(SeedsG) numel(PeriphG) numel(annotations) numel(indout)];
        SelDiseases.periph_genes{dd,1} = PeriphG;
    end
    clear SeedsG PeriphG indout check annotations Matrix gnamesD
end


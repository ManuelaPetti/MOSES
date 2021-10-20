function optimized_clustering(DB)

% Input:
% DB: databases, e.g. GO, KEGG

matrdir =   sprintf('%s Matrices WSs-CSs';
savedir =   'Results MOSES';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

load('Diseases_WSsCSs.mat')

for i=1:numel(SelDiseases.ID)
    if exist(fullfile(matrdir,sprintf('%s_%s_Matrix.mat',SelDiseases.ID{i},DB)))
        %% loading data
        GeneMatrix = loadname(fullfile(matrdir,sprintf('%s_%s_Matrix.mat',SelDiseases.ID{i},DB)));
        
        Matrix = GeneMatrix.Matrix;
        check = sum(Matrix,2);
        
        Matrix(find(check==0),:) = [];
        geneIDs = GeneMatrix.geneIDs;
        geneIDs(find(check==0),:) = [];
        SeedInAnnot = intersect(geneIDs,GeneMatrix.diseasegenes);
        SeedInAnnotPerc(i,1) = numel(SeedInAnnot)/numel(GeneMatrix.diseasegenes);
        SeedInAnnotNum(i,1) = numel(SeedInAnnot);
        
        k = 1;
        idx = ones(size(Matrix,1),1);
        Clst_SPerc = 1;
        while max(Clst_SPerc)>=.9
            k = k+1;
            idx = kmeans(Matrix,k,'Distance','hamming','Replicates',50);
            for cc=1:numel(unique(idx))
                Clst{cc,:} = geneIDs(find(idx==cc));
                ClstNodesN(1,cc) = numel(Clst{cc,:});
                Clst_SPerc(1,cc) = length(intersect(Clst{cc,:},SeedInAnnot))/numel(SeedInAnnot);
            end
        end
        MOSES.DiseaseName{i,1} = SelDiseases.Name{i,1};
        MOSES.DiseaseID{i,1} = SelDiseases.ID{i,1};
        MOSES.SelectedK (i,1) = k;
        [MaxV, indMV] = max(Clst_SPerc);
        MOSES.ClsSeedPerc(i,1) = MaxV;
        MOSES.ClustSeedID{i,:} = intersect(Clst{cc,:},SeedInAnnot);
        MOSES.PutGenesN(i,1) = ClstNodesN(1,indMV);
        MOSES.PutGenesID{i,:} = Clst{indMV,:};
        clear ClstNodesN Clst_SPerc Clst MaxV indMV
    end
end
save(fullfile(savedir,sprintf('MOSES_%s.mat',DB)),'MOSES')

 
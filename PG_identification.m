function PG_identification(DBs)

% Input:
% DBs: databases list, e.g. {'GO' 'KEGG'}

load 'Diseases_WSsCSs.mat'
savedir = 'Putative genes';
if ~exist(savedir,'dir')
    mkdir(savedir);
end

for ii=1:numel(MOSES.PutGenesID)
    clear PutGenes DB_PutGenes runIntersect
    for jj=1:length(DBs)
        load(strcat('Results MOSES\MOSES_',DBs{jj},'.mat'))
        
        if MOSES.ClsSeedPerc(ii,1)>.6
            DB_PutGenes{jj} = MOSES.PutGenesID{ii,1};
        end
    end
    runIntersect = DB_PutGenes{1};
    for kk = 2:length(PutGenes)
        runIntersect = intersect(runIntersect,DB_PutGenes{kk});
    end
    
    seeds = SelDiseases.seed{ii};
    PutGenes = setdiff(runIntersect,seeds);
    
    fid = fopen(fullfile('\Putative genes',...
        sprintf('%s.txt',SelDiseases.ID{ii,1})),'w');
    fprintf(fid,'%d\n',PutGenes);
    fclose(fid);
    
end
end

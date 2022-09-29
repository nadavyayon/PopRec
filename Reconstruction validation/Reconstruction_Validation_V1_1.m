%Reconstruction Estimation  -  19.1.21 - Nadav Yayon  
function Reconstruction_Validation_V1_1(Path,CellsFromEachSample,Gague,Vaa3DPath)
% Path = 'E:\PROJECTS\EMCU\Light_Sheet\Nadav\whiskerDeprivation4th\MouseA5\PopRec\SomaToCoubeOutput'
% init
Gague.Value = 0;
analysisFolder = 'F:\My Drive\Articles\Article_Reconstruction\Manuscript\submission_Plos_Bio\Codes';
addpath(genpath(analysisFolder))
clc 
close all
addpath(genpath(Path))
% Vaa3DDir = ' '; % path to Vaa3D exe
quest = 'Does This Neuron Make Sense?';
% SamplesPath = ' '; %path to SWC
% will go to each folder to get information on all swc files 
Global = 1; % 1 - if global coordinates, 0 - local coordinates 
zScale = 1; SomaDiam = 8; DendDiam = 1;
Directories = dir(Path);
DirDirectories = find([Directories.isdir]);
DirDirectories = DirDirectories(3:end);
counter = 0;
for i = DirDirectories
    counter = counter+1;
    Samples{counter} = Directories(i).name;
end

%% Estimate Reconstructions
cd(Path)
% see if AllCellVerdict exists on part of the data 
try         
    load('AllcellVerdict.mat');
catch 
     AllcellVerdict = cell(CellsFromEachSample*length(Samples),3);
end
counter = length(cell2mat(AllcellVerdict(:,2)))+1;

S = randperm(length(Samples)); % randomize samples 
for s = S
    cd([Path,filesep,Samples{s}]);
    load COMCoordinates
    clear Cells
    for i=1:size(COMCoordinates,1)
        Cells{i,1} = ['Cell_',num2str(i,'%03.f')];
    end
   
    cd([Path,filesep,Samples{s},'\SWC\Sorted']);
    ReconstructionConditions = length(dir)-2;
    Folder = cd;
    FileList = dir(fullfile(Folder, '**', '*.swc'));
    SizeLim = [20000 150000]; % set boundaries of expected cell sizes 
    % set file size filter and copt files to a new directory
    SelectedCells = randperm(length(Cells));
    i = 1; Flag = 0;iii=0;
    while (iii < CellsFromEachSample) && (Flag==0) && (counter<=CellsFromEachSample*length(Samples))
       i = i+1;
       
       try 
        CellIndexs = find(contains({FileList.name}',Cells{SelectedCells(i)}))
       catch me 
       end 
        CellBytes = [FileList(CellIndexs).bytes]';
       AboveBytes = (CellBytes>SizeLim(1)).*(CellBytes<SizeLim(2));
       if sum(AboveBytes)>(ReconstructionConditions/2)
           iii = iii+1;
           f = find(AboveBytes)';
           ii = randperm(length(f),1);
           
           SWCPathTemp = [FileList(CellIndexs(f(ii))).folder,filesep,FileList(CellIndexs(f(ii))).name];
           system([Vaa3DPath,'\vaa3d_msvc.exe /i ',SWCPathTemp,' &']);
           system([Vaa3DPath,'\vaa3d_msvc.exe /i ',Path,filesep,Samples{s},filesep,Cells{SelectedCells(i)},'.tif &']);
           pause(1.5);
           answer = TreeQuestion;
           answer = convertCharsToStrings(answer);
           AllcellVerdict{counter,1} = SWCPathTemp;
           switch answer
               case 'Yes' 
                   AllcellVerdict{counter,2} = 1;
               case 'No'
                   AllcellVerdict{counter,2} = 0;
               case 'Abort'
                   Flag = 1;
           end
           
            system('taskkill /im vaa3d_msvc.exe /t');
            system('taskkill /im cmd.exe /t')   
            system('taskkill /im vaa3d_msvc.exe /t');
            system('taskkill /im cmd.exe /t')   
 
            save([Path,'/','AllcellVerdict.mat'],'AllcellVerdict');
            clc

           counter = counter + 1;
           Gague.Value = 100*(counter/CellsFromEachSample/length(S));
       end    
    end
    if Flag; break; end
    if Gague.Value >= 100; break; end
end
    
    
    


%% save training files to new directory with extended name that includes origin 
Seps = strfind(Path,filesep);
PrevFolder = Path(1:Seps(end));
mkdir([PrevFolder,'TrainingSet'])
for i=1:length(AllcellVerdict)
    SepInd = strfind(AllcellVerdict{i,1},filesep);
    SampleRange = [SepInd(end-4)+1,SepInd(end-3)-1];
    CellInd = strfind(AllcellVerdict{i,1},'Cell_');
    SWCName = AllcellVerdict{i,1}(CellInd:end);
    SampleName = AllcellVerdict{i,1}(SampleRange(1):SampleRange(2));
    copyfile(AllcellVerdict{i,1},[PrevFolder,'TrainingSet\',SampleName,'_',SWCName]);
    AllcellVerdict{i,3} = [SampleName,'_',SWCName];
    AllcellVerdict{i,4} = [Path,filesep,SampleName,filesep,AllcellVerdict{i,1}(CellInd:CellInd+7),'.tif'];

end
save([PrevFolder,'AllcellVerdict.mat'],'AllcellVerdict');

%% self test for consistency 

cd(PrevFolder);
counter = 0;
Gague.Value = 100*(counter/CellsFromEachSample/length(Samples));

load AllcellVerdict.mat
if size(AllcellVerdict,2) ~= 4 
    CurrentCell = length(cell2mat(AllcellVerdict(:,5)))+1;
else
    CurrentCell = 1;
end
    counter = CurrentCell;
for Cell = CurrentCell:length(AllcellVerdict)  
    
   system([Vaa3DPath,filesep,'vaa3d_msvc.exe /i ',AllcellVerdict{Cell,1},' &']);
   system([Vaa3DPath,filesep,'vaa3d_msvc.exe /i ',AllcellVerdict{Cell,4},' &']);
   pause(2.5);
   answer = TreeQuestion;
   answer = convertCharsToStrings(answer);
   switch answer
       case 'Yes' 
           AllcellVerdict{Cell,5} = 1;
       case 'No'
           AllcellVerdict{Cell,5} = 0;
       case 'Abort'
           break
   end
   system('taskkill /im vaa3d_msvc.exe /t');
   system('taskkill /im cmd.exe /t')            
   save([PrevFolder,'/','AllcellVerdict.mat'],'AllcellVerdict');
   clc
   counter = counter + 1;
   Gague.Value = 100*(counter/CellsFromEachSample/length(Samples));
end
%% Sum Rounds and remove duplicates 
AllcellVerdictCat = [AllcellVerdict(1:end,1), AllcellVerdict(1:end,3), num2cell(cell2mat(AllcellVerdict(:,2))+cell2mat(AllcellVerdict(:,5)))];
[~ ,UInd] = unique(string(AllcellVerdictCat(:,2)),'rows');
AllcellVerdictCat = AllcellVerdictCat(UInd,:);
save([PrevFolder,'AllcellVerdictCat.mat'],'AllcellVerdictCat');
 

%% FIx3D and prune inter-nodes of Training set

cd([PrevFolder,'TrainingSet'])
DirInfo = dir('*.swc');
mkdir('N3DFixed');
parfor ii = 1:length(DirInfo)
   txt = string([Vaa3DPath,filesep,'vaa3d_msvc.exe /x N3DFix /f N3DFix /i ',DirInfo(ii).folder,filesep,DirInfo(ii).name]);
   system(txt); 
   movefile([DirInfo(ii).name,'_N3DFix.swc'],'N3DFixed');
   delete([DirInfo(ii).name,'_N3DFix_report_25_10.txt']);
end
cd('N3DFixed');
DirInfo = dir('*.swc');
mkdir('Pruned');
parfor iii = 1:length(DirInfo)
    txt = string([Vaa3DPath,filesep,'vaa3d_msvc.exe /x inter_node_pruning /f pruning /i ',DirInfo(iii).folder,filesep,DirInfo(iii).name]);
    system(txt); 
    movefile([DirInfo(iii).name,'_pruned.swc'],'Pruned');
end



%% Get Features for training 
StatPath = ['E:\PROJECTS\EMCU\CompiledNeuroM'];
text1 = [StatPath,filesep,'morph_stats1.exe ',PrevFolder,filesep,'TrainingSet',filesep,'N3DFixed',filesep,'Pruned -C ',StatPath,filesep,'config_all_13_02_19.txt -o ',PrevFolder,filesep,'FeaturesN3DFixedPruned.csv'];
system(text1);
text2 = [StatPath,filesep,'morph_stats1.exe ',PrevFolder,filesep,'TrainingSet',filesep,'N3DFixed',filesep,'Pruned -C ',StatPath,filesep,'Config_Trunk_Scholl.txt -o ',PrevFolder,filesep,'TrunkShollN3DFixedPruned.csv'];
system(text2);

%% generate data for normalization with regression learner
cd([PrevFolder])
load AllcellVerdictCat.mat
load trainedRegressionModel.mat
FeaturesTraining = readtable('FeaturesN3DFixedPruned.csv');
FeaturesTrunkSholl = importfileTrunkSholl('TrunkShollN3DFixedPruned.csv');
FeaturesTraining = [FeaturesTraining FeaturesTrunkSholl(:,[2 3 5])];
FeaturesTraining = sortrows(FeaturesTraining,'name','ascend');
[~, ind] = sort(AllcellVerdictCat(:,2));
TempTable = table('Size',[length(ind),2],'VariableTypes',{'string','double'},'VariableNames',{'NameVerdict','Rank'})
TempTable(:,1) = AllcellVerdictCat(ind,2); TempTable(:,2) = AllcellVerdictCat(ind,3);
FeaturesTrainingCategory = [FeaturesTraining TempTable];
FeaturesTrainingCategoryNorm = normalize(FeaturesTrainingCategory(:,[2:63 end]));
%% after importing, normalize data 
% sampleSize = size((FeaturesTrainingCategory),1);
% xmean = mean(FeaturesTrainingCategory{:,2:end});
% xstd = std(FeaturesTrainingCategory{:,2:end});
% [~, SortOrder] = sort(FeaturesTrainingCategory.Category); 
% Tab = tabulate(FeaturesTrainingCategory.Category) % get proportions 
% MinGroupSize = min(Tab(:,2));
% FeaturesTrainingCategory = FeaturesTrainingCategory(SortOrder,:);
% 
% % replace NaN
% FeaturesTrainingCategoryNoNan = FeaturesTrainingCategory;
% for i = 2:width(FeaturesTrainingCategory)-1
%     Nans = (isnan(FeaturesTrainingCategory{:,i}));
%     if ~isempty(find(Nans))
%         FeaturesTrainingCategoryNoNan{Nans,i} = mean(FeaturesTrainingCategory{~Nans,i});
%     end
% end
% % FeaturesTrainingCategoryStandard = (FeaturesTrainingCategory{:,2:end} - xmean(ones(sampleSize,1),:))./xstd(ones(sampleSize,1),:);
% FeaturesTrainingCategoryStandard = FeaturesTrainingCategoryNoNan;
% FeaturesTrainingCategoryStandard{:,2:end-1} = zscore(FeaturesTrainingCategoryNoNan{:,2:end-1});
% 
% Selected2s = randperm(sum(FeaturesTrainingCategory.Category==2),MinGroupSize);
% Selected1s = randperm(sum(FeaturesTrainingCategory.Category==1),MinGroupSize);
% Selected0s = randperm(sum(FeaturesTrainingCategory.Category==0),MinGroupSize);
% 
% Selected2S = FeaturesTrainingCategoryStandard(FeaturesTrainingCategory.Category==2,:);
% Selected1S = FeaturesTrainingCategoryStandard(FeaturesTrainingCategory.Category==1,:);
% Selected0S = FeaturesTrainingCategoryStandard(FeaturesTrainingCategory.Category==0,:);
% 
% Selected2SMin = Selected2S(Selected2s,:);
% Selected1SMin = Selected1S(Selected1s,:);
% Selected0SMin = Selected0S(Selected0s,:);
% 
% SelectedFetures = [2:21 23 25:48 50 51 53:58];
% TrainigData = [Selected0SMin; Selected1SMin; Selected2SMin];
% TrainigData = TrainigData(:,[SelectedFetures 62]);
% 
% 
% yfit = trainedModelLinearRegression.predictFcn(FeaturesTrainingCategoryStandard(:,SelectedFetures));
% yCat = CategorialModel.predictFcn(FeaturesTrainingCategoryStandard(:,SelectedFetures));
% Tab0 = tabulate(yCat(FeaturesTrainingCategory.Category==0)); 
% Tab1 = tabulate(yCat(FeaturesTrainingCategory.Category==1)); 
% Tab2 = tabulate(yCat(FeaturesTrainingCategory.Category==2)); 
% 
% plot(1:length(yfit),yfit,'.');

%% 3DFix and Prune inter nodes for all data
for s = 1:length(Samples)
    cd([Path,filesep,Samples{s}])
    sortedDir{s} = [Path,filesep,Samples{s},'\SWC\Sorted'];
    cd(sortedDir{s})
    FoldersPlus = dir;
    Folders{s} = FoldersPlus(logical([0 0 FoldersPlus(3:end).isdir]));
    subFolders = Folders{s};
    for i = 1:length(subFolders)
        cd(subFolders(i).name);
        DirInfo = dir('*.swc');
        mkdir('N3DFixed');
        parfor ii = 1:length(DirInfo)
           txt = string([Vaa3DPath,filesep,'vaa3d_msvc.exe /x N3DFix /f N3DFix /i ',DirInfo(ii).folder,filesep,DirInfo(ii).name]);
           system(txt); 
           movefile([DirInfo(ii).name,'_N3DFix.swc'],'N3DFixed');
           delete([DirInfo(ii).name,'_N3DFix_report_25_10.txt']);
        end
        cd('N3DFixed');
        DirInfo = dir('*.swc');
        mkdir('Pruned');
        parfor iii = 1:length(DirInfo)
            txt = string([Vaa3DPath,filesep,'vaa3d_msvc.exe /x inter_node_pruning /f pruning /i ',DirInfo(iii).folder,filesep,DirInfo(iii).name]);
            system(txt); 
            movefile([DirInfo(iii).name,'_pruned.swc'],'Pruned');
        end
        cd(sortedDir{s})
    end
    
end

%% FIx3D and prune inter-nodes of all cells (this might take some time)

% import boostedtree
% get all adresses and features 
for s = 1:length(Samples)
    cd([Path,'/',Samples{s}])
    Cells = struct2cell(dir('*.tif'))';
    CellNames{s} = Cells(:,1);
    sortedDir{s} = [Path,filesep,Samples{s},'\SWC\Sorted'];
    cd(sortedDir{s})
    mkdir('RelevantCells');
    FoldersPlus = dir;
    Folders{s} = FoldersPlus(logical([0 0 FoldersPlus(3:end).isdir])); 
    AllNames = [];
    AllCellSize = zeros(length(CellNames{s}),length(Folders));
    subFolders = Folders{s};
    mkdir('RelevantCells');
    for i = 1:length(subFolders)-1
        eval([subFolders(i).name,'_SWCInfo = dir(''',subFolders(i).name,'','\N3DFixed\Pruned'')'])
        temp = struct2cell(eval([subFolders(i).name,'_SWCInfo']))';
        IndexXYZ = find(contains(temp(:,1),'XYZ'));
        IndexAll = ones(length(temp),1);
        IndexNotXYZ = IndexAll;
        IndexNotXYZ([1; 2; IndexXYZ]) = 0; 
        IndexNotXYZ = find(IndexNotXYZ);
        if Global 
           Index = IndexXYZ;
        else 
           Index = IndexNotXYZ;
        end
        relevantIndex =  find((cell2mat(temp(Index(:),4)) > SizeLim(1)) .*  (cell2mat(temp(Index(:),4)) < SizeLim(2)));
        relevantCells = Index(relevantIndex);
        for Cell = 1:length(relevantCells)
         copyfile([subFolders(i).name,'\N3DFixed\Pruned\',temp{relevantCells(Cell),1}],[sortedDir{s},'\RelevantCells']);       
        end
    end
    formatOut = 'mmddyyyy';
    Date = datestr(now,formatOut);
    text1 = [StatPath,filesep,'morph_stats1.exe ',sortedDir{s},filesep,'RelevantCells -C ',StatPath,filesep,'config_all_13_02_19.txt -o ',sortedDir{s},filesep,'FeaturesN3DFixedPruned_',Date,'.csv'];
    system(text1);
    text2 = [StatPath,filesep,'morph_stats1.exe ',sortedDir{s},filesep,'RelevantCells -C ',StatPath,filesep,'Config_Trunk_Scholl.txt -o ',sortedDir{s},filesep,'FeaturesTrunkShollN3DFixedPruned_',Date,'.csv'];
    system(text2);
%     Features{s} = readtable(['Features_',Date,'.csv']);
%     for colomb = 2:size(Features{s},2) % replace NaNs with variable mean 
%         Features{s}{isnan(Features{s}{:,colomb}),colomb} = mean(Features{s}{find(~isnan(Features{s}{:,colomb})),colomb});
%     end
end

 %% Predict for all cells 
load([PrevFolder,filesep,'trainedRegressionModel.mat']);
clear FeaturesTrainingCategoryNorm
for s = 1:length(Samples)
    FeaturesTraining = readtable([Path,filesep,Samples{s},filesep,'SWC',filesep,'Sorted',filesep,'FeaturesN3DFixedPruned_04222021.csv']);
    FeaturesTrunkSholl = importfileTrunkSholl([Path,filesep,Samples{s},filesep,'SWC',filesep,'Sorted',filesep,'FeaturesTrunkShollN3DFixedPruned_04222021.csv']);
    FeaturesTraining = [FeaturesTraining FeaturesTrunkSholl(:,[2 3 5])];
    FeaturesTraining = sortrows(FeaturesTraining,'name','ascend');  
    FeaturesTrainingAll{s} = FeaturesTraining;
    FeaturesTrainingCategoryNorm{s} = normalize(FeaturesTraining(:,[2:63]));
end
%% Predict
figure(1)
suptitle('Prediction Ranks') 
for s = 1:length(Samples)
    subplot(6,4,s)
    title(['Sample No ',num2str(s)]);
    Prediction{s} = trainedRegressionModel.predictFcn(FeaturesTrainingCategoryNorm{s}(:,[1:59 62 60 61]));
    plot(1:1:length(Prediction{s}),Prediction{s}(1:1:end),'.');
    PredictionsT = table(Prediction{s});
    PredictionsT.Properties.VariableNames{1} = 'Predictions';
    FeaturesTrainingAllP{s} = [FeaturesTrainingAll{s} PredictionsT];
end

%% save combined features and predictions
cd(PrevFolder)
save('FeaturesTrainingAllP.mat','FeaturesTrainingAllP');


%% load Features and seperate the best cells
cd(PrevFolder)
load('FeaturesTrainingAllP.mat');


for s = 1:length(FeaturesTrainingCategoryNorm)
  
    sortedDir{s} = [Path,filesep,Samples{s},filesep,'SWC',filesep,'Sorted'];
    cd(sortedDir{s})
    cd('RelevantCells');
    try
        rmdir GoodCells s
    end
    mkdir('GoodCells');
    bestCells{s} = FeaturesTrainingAllP{s};
    bestCells{s}{:,1} = {''};
    bestCells{s}{:,2:end} = NaN;
    bestCells{s} = bestCells{s}(1,:);
    counter = 1;
    TempFeatures = FeaturesTrainingAllP{s};
    relevantCells =  find(TempFeatures{:,'Predictions'} >= 0);
    TempFeatures  = TempFeatures(relevantCells,:);
    for cells = 1:height(TempFeatures)
        % find other good cells
        if ~isempty(TempFeatures) 
            CellNameTemp = char(strjoin([TempFeatures{1,'name'},'.swc']));
            allcellswiththatname = find(not(cellfun('isempty',strfind(TempFeatures{:,'name'},CellNameTemp(1:8))))) ; 
            
            if length(allcellswiththatname)>1
                clear TempPredictions
                for i = 1:length(allcellswiththatname)
%                     str = TempFeatures{allcellswiththatname(i),'name'}{1};
%                     tempIdx = find(str == '_');
                    TempPredictions(i) = TempFeatures{allcellswiththatname(i),'Predictions'};
                  
                end
                [a] = max(TempPredictions);
                FindMax = find(TempPredictions == a);
%                 if length(FindMax) > 1
%                     [~ ,FindMinFT] = max(FT(FindMax));
%                     FindMax = FindMax(FindMax);
%                 end
            else
                FindMax = 1;
            end
                
            CellNameTemp = char(strjoin([TempFeatures{allcellswiththatname(FindMax(1)),'name'},'.swc']));
            copyfile(CellNameTemp(~isspace(CellNameTemp)),[sortedDir{s},'\RelevantCells\GoodCells']);       
            FeaturesPassThresholdCellsAll{s}(counter,:) = TempFeatures(FindMax(1),:);
            TempFeatures(allcellswiththatname,:) = [];
            counter = counter+1 ;
        end
    end
end
cd(PrevFolder)
save('FeaturesPassThresholdCellsAll.mat','FeaturesPassThresholdCellsAll')
%% fix features for sholl for best cells only and fix NaNs
cd(PrevFolder)
load('FeaturesPassThresholdCellsAll.mat')
for s = 1:length(Samples)  
    SchollTable = FeaturesPassThresholdCellsAll{s}.sholl_frequency;
    clear sholl
    for i = 1:length(SchollTable)
        clear C D E 
        C = strsplit(SchollTable(i),'[');
        D = str2double(strsplit(C(2),','));
        E = D(2:2:30);
        sholl{i} = E;
    end
    SchollTable = [table(SchollTable) table(sholl')];
    SchollTable.Properties.VariableNames(2) = {'sholl'};
    SchollTable.Properties.VariableNames(1) = {'sholl_frequency'};

    FeaturesPassThresholdCellsAllwSholl{s} = [FeaturesPassThresholdCellsAll{s} SchollTable(:,2)];
    
    for colomb = [2:(size(FeaturesPassThresholdCellsAllwSholl{s},2)-3)] % replace NaNs with variable mean 
        FeaturesPassThresholdCellsAllwSholl{s}{isnan(FeaturesPassThresholdCellsAllwSholl{s}{:,colomb}),colomb} = mean(FeaturesPassThresholdCellsAllwSholl{s}{find(~isnan(FeaturesPassThresholdCellsAllwSholl{s}{:,colomb})),colomb});
        isnan(FeaturesPassThresholdCellsAllwSholl{s}{:,colomb})
    end
end
save('FeaturesPassThresholdCellsAllwSholl.mat','FeaturesPassThresholdCellsAllwSholl')

%% add soma coordinates and write to  new tabel

for s = 1:length(Samples)
       Coordinates = table([0],[0],[0],[0],'VariableNames',{'Soma_x','Soma_y','Soma_z','SomaDiameter'});
       cd([Path,filesep,Samples{s}])
       load COMCoordinates.mat;
       load SomaStats.mat;
       

        for i = 1:length(FeaturesPassThresholdCellsAllwSholl{s}{:,'name'})
            CellName = [FeaturesPassThresholdCellsAllwSholl{s}{i,'name'}{1},'.swc'];
            CellNumber = str2num(CellName(6:8));
            Coordinates{i,:} = [COMCoordinates(CellNumber,:) SomaStats(CellNumber,:).EquivDiameter] ;
        end
        cd(PrevFolder)
        writetable([FeaturesPassThresholdCellsAllwSholl{s} Coordinates],['AllFeaturesSelectedCells_',Samples{s},'.csv'])
        
 end 
end

%%%%%%%%%%%%%%%%%%%%%%



function TrunkShollN3DFixedPruned = importfileTrunkSholl(filename, dataLines)
%IMPORTFILE Import data from a text file
%  TRUNKSHOLLN3DFIXEDPRUNED = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as a table.
%
%  TRUNKSHOLLN3DFIXEDPRUNED = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  TrunkShollN3DFixedPruned = importfile("E:\PROJECTS\EMCU\4thWhiskerDeprived\SomaToCube\TrunkShollN3DFixedPruned.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 20-Apr-2021 22:09:04

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["name", "allmean_trunk_angle", "allstd_trunk_angle", "mean_soma_radius", "sholl_frequency"];
opts.VariableTypes = ["string", "double", "double", "double", "string"];

% Specify file level properties
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["name", "sholl_frequency"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["name", "sholl_frequency"], "EmptyFieldRule", "auto");

% Import the data
TrunkShollN3DFixedPruned = readtable(filename, opts);

end
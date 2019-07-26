%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_ESCLogicSystemAnalysis.m
%%%%
%%%%        - Required Info:
%%%%            - DPA Boundary Information
%%%%            - DPA Inner Boundary Information
%%%%            - DPA Grid Point Information
%%%%
%%%%        - Function that Generates and Saves ESC Logic System Info for the given DPA
%%%%            - ESC sensor Information
%%%%            - DPA Grid Point Detection Information
%%%%            - DPA Grid Point Coverage Information
%%%%            - ESC sensor Combinations Priority List; along with
%%%%                    - Principal Sensor Combination Prioritu List
%%%%                    - Neighborhing Sensor Combination Priority List
%%%%                    - Cross Sensor Combination Priority List
%%%%            - ESC sensor DPA Coverage Combinations Priority List
%%%%                    - along with Successful Sensor Combination Information
%%%%
%%%%        - Generates and saves DPA Coverage Combinations Table
%%%%            - ESCLogicSystem_CombinationSensorSiteID.xlsx
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% DPA Information upto 12-13-2017
% east_dpa_1.geojson (1 to 26);
% east_dpa_bath.geojson (27)
% east_dpa_mayport.geojson (28)
% east_dpa_norfolk.geojson (29)
% east_dpa_pascagoula.geojson (30)
%
% west_dpa_1.geojson (1 to 14);
% west_dpa_Alameda.geojson (15)
% west_dpa_bremerton.geojson (16)
% west_dpa_LongBeach.geojson (17)
% west_dpa_SanDiego.geojson (18)
% west_dpa_Washington.geojson (19)
%
% alaska_dpa_1.geojson (1 to 38)
% hawaii_dpa_1.geojson (1 to 9)
% hawaii_dpa_PearlHarbor.geojson (10)
% puerto_rico_ver2_dpa_1.geojson (1 to 5)
% Guam_dpa_1.geojson (1 and 2)
% Guam_dpa_ApraHarbor.geojson (3)






function Function_ELSAnalysis_ST2(DPA_Zone, DPA_ID)



format long;




if strcmp(DPA_Zone, 'east') && DPA_ID == 27
    DPA_IDstr = 'Norfolk';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 28
    DPA_IDstr = 'Mayport';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 29
    DPA_IDstr = 'Pascagoula';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 30
    DPA_IDstr = 'Pensacola';
elseif strcmp(DPA_Zone, 'east') && DPA_ID == 31
    DPA_IDstr = 'WebsterField';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 15
    DPA_IDstr = 'SanDiego';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 16
    DPA_IDstr = 'Alameda';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 17
    DPA_IDstr = 'Long-Beach';
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 18
    DPA_IDstr = 'Bremerton-Everett';
else
    DPA_IDstr = DPA_ID;
end     % If ends after assigning numeric values to DPA with string names






% Generate the DPA Name
DPA_Name = strcat(DPA_Zone, '_', num2str(DPA_IDstr));
DPAfolderName = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));



ELS_dataSaveFolder = strcat('C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Data_20190531\', DPAfolderName);
if exist(ELS_dataSaveFolder, 'dir') ~= 7
    mkdir(ELS_dataSaveFolder)
end



DPArelated_dataSaveFolder = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Original DPA and ESC Data';
DPAspecific_dataSaveFolder = strcat(DPArelated_dataSaveFolder, '\', DPAfolderName);












xx=1;   % Completes Input initialization
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPA Grid Point Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DPA_gridInfo_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_gridInfo.mat');
DPA_gridInfo_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_gridInfo_fileName);
if exist(DPA_gridInfo_fullFileName, 'file') == 2
    load(DPA_gridInfo_fullFileName, 'DPA_gridInfo')
else
    DPA_gridInfo = Function_GenerateDPAGridPoints(DPA_Zone, DPA_ID);
end

gridPointCount = size(DPA_gridInfo, 1);
DPA_gridInfo(:,4) = DPA_gridInfo(:,3); 
DPA_gridInfo(:,3) = (1:gridPointCount)';



% Save DPA_gridInfo
DPA_gridInfo_fileName = strcat(DPA_Name, '_DPA_gridInfo.mat');
DPA_gridInfo_fullFileName = fullfile(ELS_dataSaveFolder, DPA_gridInfo_fileName);
if exist(DPA_gridInfo_fullFileName, 'file') == 2
    delete(DPA_gridInfo_fullFileName)
end
save(DPA_gridInfo_fullFileName, 'DPA_gridInfo')
clear DPA_gridInfo_fileName 




%%% Create the DPA Grid Point Info table to be saved as .xlsx file
GridPointIndex = (1:size(DPA_gridInfo,1))'; 
DPAGridPointInfoTable = table(GridPointIndex, DPA_gridInfo(:,1), DPA_gridInfo(:,2), DPA_gridInfo(:,4));
DPAGridPointInfoTable.Properties.VariableNames{1} = 'GridPointIndex';
DPAGridPointInfoTable.Properties.VariableNames{2} = 'GridPointLon';
DPAGridPointInfoTable.Properties.VariableNames{3} = 'GridPointLat';
DPAGridPointInfoTable.Properties.VariableNames{4} = 'InnerBoundaryIndex';

% Save the ESC sensor combination Information (.xlsx file)


DPAGridPointInfo_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_DPAGridPointInfo.xlsx');
DPAGridPointInfo_fullFileName = fullfile(ELS_dataSaveFolder, DPAGridPointInfo_filename);
if exist(DPAGridPointInfo_fullFileName, 'file') == 2
    delete(DPAGridPointInfo_fullFileName)
end

writetable(DPAGridPointInfoTable, DPAGridPointInfo_fullFileName, 'Sheet', 1, 'Range','A1')
clear DPAGridPointInfo_fullFileName








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPA Boundary Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DPA_Boundary_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_Boundary.mat');
DPA_Boundary_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_Boundary_fileName);
if exist(DPA_Boundary_fullFileName, 'file') == 2
    if strcmp(DPA_Zone, 'west') && DPA_ID == 13
        load(DPA_Boundary_fullFileName, 'DPA_Boundary_cell')
        
        % Save DPA_Boundary
        DPA_Boundary_fullFileName = fullfile(ELS_dataSaveFolder, DPA_Boundary_fileName);
        if exist(DPA_Boundary_fullFileName, 'file') == 2
            delete(DPA_Boundary_fullFileName)
        end
        save(DPA_Boundary_fullFileName, 'DPA_Boundary_cell')
        
    else
        load(DPA_Boundary_fullFileName, 'DPA_Boundary')
        
        % Save DPA_Boundary
        DPA_Boundary_fullFileName = fullfile(ELS_dataSaveFolder, DPA_Boundary_fileName);
        if exist(DPA_Boundary_fullFileName, 'file') == 2
            delete(DPA_Boundary_fullFileName)
        end
        save(DPA_Boundary_fullFileName, 'DPA_Boundary')
    end
    
else
    if strcmp(DPA_Zone, 'west') && DPA_ID == 13
        [~, DPA_Boundary_cell] = Function_GenerateDPABoundary(DPA_Zone, DPA_ID);        
        
        DPA_Boundary_fullFileName = fullfile(ELS_dataSaveFolder, DPA_Boundary_fileName);
        if exist(DPA_Boundary_fullFileName, 'file') == 2
            delete(DPA_Boundary_fullFileName)
        end
        save(DPA_Boundary_fullFileName, 'DPA_Boundary_cell')
        
    else
        DPA_Boundary = Function_GenerateDPABoundary(DPA_Zone, DPA_ID);
        
        DPA_Boundary_fullFileName = fullfile(ELS_dataSaveFolder, DPA_Boundary_fileName);
        if exist(DPA_Boundary_fullFileName, 'file') == 2
            delete(DPA_Boundary_fullFileName)
        end
        save(DPA_Boundary_fullFileName, 'DPA_Boundary')
        
    end
    
end
clear DPA_Boundary_fileName DPA_Boundary_fullFileName














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPA Inner Boundary Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DPA_innerBoundary_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPA_innerBoundary.mat');
DPA_innerBoundary_fullFileName = fullfile(DPAspecific_dataSaveFolder, DPA_innerBoundary_fileName);
if exist(DPA_innerBoundary_fullFileName, 'file') == 2
    load(DPA_innerBoundary_fullFileName, 'DPA_innerBoundary')
    
    % Save DPA_innerBoundary
    DPA_innerBoundary_fullFileName = fullfile(ELS_dataSaveFolder, DPA_innerBoundary_fileName);
    if exist(DPA_innerBoundary_fullFileName, 'file') == 2
        delete(DPA_innerBoundary_fullFileName)
    end
    save(DPA_innerBoundary_fullFileName, 'DPA_innerBoundary')    
    
else
    DPA_innerBoundary = Function_GenerateDPAInnerBoundary(DPA_Zone, DPA_ID);
    
    DPA_innerBoundary_fullFileName = fullfile(ELS_dataSaveFolder, DPA_innerBoundary_fileName);
    if exist(DPA_innerBoundary_fullFileName, 'file') == 2
        delete(DPA_innerBoundary_fullFileName)
    end
    save(DPA_innerBoundary_fullFileName, 'DPA_innerBoundary')
    
end
clear DPA_innerBoundary_fileName DPA_innerBoundary_fullFileName








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESC Sensor Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ESCsensorInfo = Function_GenerateESCsensorInformation(DPA_Zone, DPA_ID);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESC Sensor Antenna Gain Pattern Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Function_GenerateAntennaPattern_MultipleESC(DPA_Zone, DPA_ID, ESCsensorInfo)








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Go TO THE PYTHON CODE TO GENERATE THE GRID POINT DETECTION AND COVERAGE INFORMATION
%%%     Grid Point Detection Info: Return with the GridPointDetectionInfo.xlsx file
%%%     Grid Point Coverage Info: Return with the GridPointCoverageInfo.xlsx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














xx=1;       fprintf('\nESC Logic System: Grid Point Detection and Coverage Information :: Use Python code\n')      % Completes gathering ESC sensor information

pause














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DPA Grid Point Detection and Coverage Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELS_GridPointDetectionInfo_fullFileName = strcat(ELS_dataSaveFolder, '\', DPA_Zone, '_', num2str(DPA_IDstr), '_DPAGridPointDetectionInfo.xlsx');
GridPointDetectionInfo_Data = readtable(ELS_GridPointDetectionInfo_fullFileName, 'sheet', 'Sheet1');


GridPointDetectionInfo_DataArray = table2array(GridPointDetectionInfo_Data);
DetectionReliability = max(GridPointDetectionInfo_DataArray(:,2:end),[],2);
GridPointswithNoDetection = find(DetectionReliability==0); %#ok<NASGU>








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate the Grid Point Coverage set for each of the ESC sensors 
%%%     GridPointCoverage: for each of the ESC sensors (N) associated with the selected DPA
%%%         cell{for each of the Grid Points (Mi) covered by the ith ESC sensor: 
%%%         [GridPintIndex, MaxReliabilityAchieved]: (Nx1 cell -> (Mix2 Double))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nESC Logic System: Grid Point detection and coverage info generation :: Start \n')


% Some of the West DPAs still do not have required coverage
if strcmp(DPA_Zone, 'west') && DPA_ID==1
    
    Inside75km_ReliabilityThreshold = 0.75;
    Beyond75km_ReliabilityThreshold = 0.00;
    
elseif strcmp(DPA_Zone, 'west') && ismember(DPA_ID, [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
    
     Inside75km_ReliabilityThreshold = 0.95;
     Beyond75km_ReliabilityThreshold = 0.00;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==26
    
    Inside75km_ReliabilityThreshold = 0.95;
    Beyond75km_ReliabilityThreshold = 0.00;
    
else
    
    Inside75km_ReliabilityThreshold = 0.95;
    Beyond75km_ReliabilityThreshold = 0.50;
    
end




escCount = size(ESCsensorInfo.sensorIndex,1);
GridPointDetectionInfo = cell(escCount,1);
GridPointCoverageInfo = cell(escCount,1);
for escIndex = 1 : escCount     % For each ESC sensors
    
    
    % Get the Grid Point Detection and Coverage Information
    if escIndex==1
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x0;
    elseif escIndex==2
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x1;
    elseif escIndex==3
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x2;
    elseif escIndex==4
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x3;
    elseif escIndex==5
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x4;
    elseif escIndex==6
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x5;
    elseif escIndex==7
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x6;
    elseif escIndex==8
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x7;
    elseif escIndex==9
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x8;
    elseif escIndex==10
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x9;
    elseif escIndex==11
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x10;
    elseif escIndex==12
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x11;
    elseif escIndex==13
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x12;
    elseif escIndex==14
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x13;
    elseif escIndex==15
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x14;
    elseif escIndex==16
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x15;
    elseif escIndex==17
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x16;
    elseif escIndex==18
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x17;
    elseif escIndex==19
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x18;
    elseif escIndex==20
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x19;
    elseif escIndex==21
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x20;
    elseif escIndex==22
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x21;
    elseif escIndex==23
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x22;
    elseif escIndex==24
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x23;
    elseif escIndex==25
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x24;
    elseif escIndex==26
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x25;
    elseif escIndex==27
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x26;
    elseif escIndex==28
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x27;
    elseif escIndex==29
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x28;
    elseif escIndex==30
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x29;
    elseif escIndex==31
        currentESC_GridPointDetectionInfo = GridPointDetectionInfo_Data.x30;
    end
    
    
    
    detectedGPindex = 0;        coveredGPindex = 0;
    for gridpointIndex = 1 : gridPointCount      % For each DPA grid point
        
%         if currentESC_GridPointDetectionInfo(gridpointIndex,1)>0
        if currentESC_GridPointDetectionInfo(gridpointIndex,1)>=0    
            detectedGPindex = detectedGPindex + 1;
            GridPointDetectionInfo{escIndex, 1}(detectedGPindex, 1) = gridpointIndex;  % record the index of the detected grid point
            GridPointDetectionInfo{escIndex, 1}(detectedGPindex, 2) = currentESC_GridPointDetectionInfo(gridpointIndex,1);  % record the maximum detection reliability
            
        end
        
        
        if DPA_gridInfo(gridpointIndex,4)==1
            
            if ge(currentESC_GridPointDetectionInfo(gridpointIndex,1), Inside75km_ReliabilityThreshold)
                
                coveredGPindex = coveredGPindex + 1;
                GridPointCoverageInfo{escIndex, 1}(coveredGPindex, 1) = gridpointIndex;  % record the index of the covered grid point
                GridPointCoverageInfo{escIndex, 1}(coveredGPindex, 2) = currentESC_GridPointDetectionInfo(gridpointIndex,1);  % record the maximum coverage reliability
            end
            
        elseif DPA_gridInfo(gridpointIndex,4)==0
            
            if ge(currentESC_GridPointDetectionInfo(gridpointIndex,1), Beyond75km_ReliabilityThreshold)
                coveredGPindex = coveredGPindex + 1;
                GridPointCoverageInfo{escIndex, 1}(coveredGPindex, 1) = gridpointIndex;  % record the index of the covered grid point
                GridPointCoverageInfo{escIndex, 1}(coveredGPindex, 2) = currentESC_GridPointDetectionInfo(gridpointIndex,1);  % record the maximum coverage reliability
            end
            
        end
        
        
        
        
        
    end     % end of DPA grid point loop
    
    
    
    if ~isempty(GridPointCoverageInfo{escIndex, 1})
        ESCsensorInfo.CoverageReliability(escIndex,1) = min(GridPointCoverageInfo{escIndex, 1}(:,2));
    else
        ESCsensorInfo.CoverageReliability(escIndex,1) = 0;
    end
    
    
    
end     % end of ESC sensor loop






xx=1;








% Save Grid Point Detection and Coverage Performance Information
GridPointDetectionandCoverageInfo_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_GridPointDetectionandCoverageInfo.mat');
GridPointDetectionandCoverageInfo_fullFileName = fullfile(ELS_dataSaveFolder, GridPointDetectionandCoverageInfo_fileName);
if exist(GridPointDetectionandCoverageInfo_fullFileName, 'file') == 2
    delete(GridPointDetectionandCoverageInfo_fullFileName)
end
save(GridPointDetectionandCoverageInfo_fullFileName, 'GridPointDetectionInfo', 'GridPointCoverageInfo')
clear GridPointDetectionandCoverageInfo_fileName GridPointDetectionandCoverageInfo_fullFileName








% Save ESCsensorInfo.mat
ESCsensorInfo_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_ESCsensorInfo.mat');
ESCsensorInfo_fullFileName = fullfile(ELS_dataSaveFolder, ESCsensorInfo_fileName);
if exist(ESCsensorInfo_fullFileName, 'file') == 2
    delete(ESCsensorInfo_fullFileName)
end
save(ESCsensorInfo_fullFileName, 'ESCsensorInfo')
clear ESCsensorInfo_fileName ESCsensorInfo_fullFileName














xx=1;       fprintf('\nESC Logic System: Grid Point detection and coverage info generation :: Complete\n')














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate ESC Sensor Combination Priority List Ordering according to SensorCategory and CoverageReliability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nESC Logic System: ESC sensor Coverage Combination Priority List :: Start\n')    

%%% Generate the list of ordered ESC sensors with individual priority
% Generate the list of ordered Principal ESC sensors with individual priority score
PrincipalESCindex = find(ESCsensorInfo.SensorCategory==1);
[~, PrincipalESCorderedIndex] = sort(ESCsensorInfo.CoverageReliability(PrincipalESCindex), 'descend');

% Generate the list of ordered Neighboring ESC sensors with individual priority score
NeighboringESCindex = find(ESCsensorInfo.SensorCategory==0);
[~, NeighboringESCorderedIndex] = sort(ESCsensorInfo.CoverageReliability(NeighboringESCindex), 'descend');
NeighboringESCorderedIndex = NeighboringESCorderedIndex + size(PrincipalESCorderedIndex,1);

% Generate the list of ordered ESC sensors with individual priority
ESCsensorOrderedIndex_temp = [PrincipalESCorderedIndex; NeighboringESCorderedIndex];
ESCsensorOrderedIndex = [ESCsensorOrderedIndex_temp (1:escCount)']; % ESCsensorOrderedIndex = [ESCindex ESCpriority]






%%% Generate All possible ESC sensor combinations and the corresponding Combination Priority Score
ESCsensorCombinations = double.empty(0, escCount+1);
for escIndex = 1 : escCount
    
    ESCsensorCombinations_temp = zeros(nchoosek(escCount, escIndex), escIndex+1);
    ESCsensorCombinations_temp(:, 2:end) = (nchoosek(ESCsensorOrderedIndex(:,1), escIndex));
    
    
    for tempIndex = 1 : nchoosek(escCount, escIndex)
        
        currentESCsensorCombination = ESCsensorCombinations_temp(tempIndex, 2:end)';        
        ESCsensorPriority1 = ESCsensorOrderedIndex(ismember(ESCsensorOrderedIndex(:,1), currentESCsensorCombination, 'rows'), 2);        
        ESCsensorCombinations_temp(tempIndex, 1) = sum(ESCsensorPriority1);        
        
    end
    
    ESCsensorCombinations_temp = [ESCsensorCombinations_temp, zeros(nchoosek(escCount,escIndex), escCount-escIndex)]; %#ok<AGROW>
    
    ESCsensorCombinations = [ESCsensorCombinations; ESCsensorCombinations_temp]; %#ok<AGROW>
    clear ESCsensorCombinations_temp
    
    
end











%%% Generate the ESC sensor combination Priority Lists
ESCsensorCombinations_withoutScore = ESCsensorCombinations(:, 2:end);

% Generate the ESC sensor Principal-Combination Priority List
PrincipalSensorCombinaitons_GlobalSet = ESCsensorCombinations(all(ismember(ESCsensorCombinations_withoutScore, [PrincipalESCindex;0]),2), :);
[~, PrincipalSensorCombinaitons_GlobalSet_PriorityIndex] = sort(PrincipalSensorCombinaitons_GlobalSet(:,1));
PrincipalSensorCombinaitons_PriorityListwithScore_GlobalSet = PrincipalSensorCombinaitons_GlobalSet(PrincipalSensorCombinaitons_GlobalSet_PriorityIndex,:);

% Generate the ESC sensor Neighboring-Combination Priority List
NeighboringSensorCombinaitons_GlobalSet = ESCsensorCombinations(all(ismember(ESCsensorCombinations_withoutScore, [NeighboringESCindex;0]),2), :);
[~, NeighboringSensorCombinaitons_GlobalSet_PriorityIndex] = sort(NeighboringSensorCombinaitons_GlobalSet(:,1));
NeighboringSensorCombinaitons_PriorityListwithScore_GlobalSet = NeighboringSensorCombinaitons_GlobalSet(NeighboringSensorCombinaitons_GlobalSet_PriorityIndex,:);

% Generate the ESC sensor Cross-Combination Priority List
CrossSensorCombinaitons_PriorityListwithScore_GlobalSetTemp = ...
    ESCsensorCombinations(~ismember(ESCsensorCombinations, PrincipalSensorCombinaitons_PriorityListwithScore_GlobalSet, 'rows'), :);
CrossSensorCombinaitons_PriorityListwithScore_GlobalSet = CrossSensorCombinaitons_PriorityListwithScore_GlobalSetTemp...
    (~ismember(CrossSensorCombinaitons_PriorityListwithScore_GlobalSetTemp, NeighboringSensorCombinaitons_PriorityListwithScore_GlobalSet, 'rows'), :);

% Generate the ESC sensor combination Priority List 
ESCsensorCombination_PriorityListwithScore_GlobalSet = [PrincipalSensorCombinaitons_PriorityListwithScore_GlobalSet; ...
    CrossSensorCombinaitons_PriorityListwithScore_GlobalSet; NeighboringSensorCombinaitons_PriorityListwithScore_GlobalSet];


% Save ESCsensorCombinationsPriorityList
ESCsensorCombination_PriorityList_GlobalSet_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_ESCsensorCombination_PriorityListwithScore_GlobalSet.mat');
ESCsensorCombination_PriorityList_GlobalSet_fullFileName = fullfile(ELS_dataSaveFolder, ESCsensorCombination_PriorityList_GlobalSet_fileName);
save(ESCsensorCombination_PriorityList_GlobalSet_fullFileName, 'ESCsensorCombination_PriorityListwithScore_GlobalSet', ...
    'PrincipalSensorCombinaitons_PriorityListwithScore_GlobalSet', 'NeighboringSensorCombinaitons_PriorityListwithScore_GlobalSet', 'CrossSensorCombinaitons_PriorityListwithScore_GlobalSet');
clear ESCsensorCombination_PriorityList_GlobalSet_fileName ESCsensorCombination_PriorityList_GlobalSet_fullFileName



















xx=1;       fprintf('\nESC Logic System: All Possible ESC Sensor Combination Priority List Generation :: Complete\n')


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate ESC Sensor Combination Priority List for the sensor combinations that satisfies DPA Coverage requirements for all Grid Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DPACoverageGridPointRequirements = DPA_gridInfo(:,3);
DPACoverageGridPointRequirements(:,2) = Beyond75km_ReliabilityThreshold * ones(gridPointCount,1);
DPACoverageGridPointRequirements(DPA_gridInfo(:,4)==1, 2) = Inside75km_ReliabilityThreshold * ones(sum(DPA_gridInfo(:,4)==1),1);

% East-26 rear grid points are in Mexico and NTIA agreed to relax requirement for those grid points
if strcmp(DPA_Zone, 'east') && (DPA_ID==26)    
    East26_rearGridPointIndex = [5, 13, 25, 41, 61]';    
    DPACoverageGridPointRequirements(East26_rearGridPointIndex,2) = 0;
end



successfulCombinationCount = 0;
SuccessfulSensorCombinaitons = zeros(size(ESCsensorCombination_PriorityListwithScore_GlobalSet(:, 2:end)));
for sensorCombinationIndex = 1 : size(ESCsensorCombination_PriorityListwithScore_GlobalSet,1)
    
    newSuccessfulCombinationFlag = 0;
    
    iterationCombination = ESCsensorCombination_PriorityListwithScore_GlobalSet(sensorCombinationIndex, 2:end)';
    iterationCombination(iterationCombination==0) = [];
    
    
    CombinationSensorCoverageReliability_All = zeros(size(DPA_gridInfo,1), size(iterationCombination,1));
    for sensorIndex = 1 : size(iterationCombination,1)
        
        if ~isempty(GridPointCoverageInfo{iterationCombination(sensorIndex)})
            CombinationSensorCoverageReliability_All(GridPointCoverageInfo{iterationCombination(sensorIndex)}(:,1), sensorIndex) ...
                = GridPointCoverageInfo{iterationCombination(sensorIndex)}(:,2);
        end
        
    end
    
    
    
    CombinationSensorCoverageReliability = max(CombinationSensorCoverageReliability_All, [], 2);
    if any(DPACoverageGridPointRequirements(:,2) > CombinationSensorCoverageReliability)
        continue        % Consider only those sensor combinations that satisfies DPA Coverage requirements
        
    else
        
        if successfulCombinationCount~=0        % Exclude all the super sets of ESC sensor combinations with higher priority
            
            for prevSuccessfulCombinationIndex = 1 : successfulCombinationCount
                
                iterationPrevSuccessfulCombination = DPACoverageCombinationPriorityList.SensorCombination{prevSuccessfulCombinationIndex,1};                
                
                if all(ismember(iterationPrevSuccessfulCombination, iterationCombination))
                    newSuccessfulCombinationFlag = 0;
                    break
                else
                    newSuccessfulCombinationFlag = 1;
                end
                
            end
            
        end
        
        
        if successfulCombinationCount==0 || newSuccessfulCombinationFlag == 1
            
            successfulCombinationCount = successfulCombinationCount + 1;
            SuccessfulSensorCombinaitons(successfulCombinationCount, 1:size(iterationCombination,1)) = iterationCombination';
            
            DPACoverageCombinationPriorityList.SensorCombination{successfulCombinationCount,1} = iterationCombination;
            DPACoverageCombinationPriorityList.CombinationSensorSiteID{successfulCombinationCount,1} = ESCsensorInfo.SiteID(iterationCombination);
            DPACoverageCombinationPriorityList.GPcoverageInfo{successfulCombinationCount,1} = [DPACoverageGridPointRequirements(:,1), CombinationSensorCoverageReliability];
            
        end
        
        
    end
    
    
    
    
    
    
    
    
end


SuccessfulSensorCombinaitons(SuccessfulSensorCombinaitons(:,1)==0, :) = [];
numberofSensorsinCombinations = zeros(size(SuccessfulSensorCombinaitons,1),1);
for combinationIndex = 1 : size(SuccessfulSensorCombinaitons,1)
    numberofSensorsinCombinations(combinationIndex,1) = find(SuccessfulSensorCombinaitons(combinationIndex,:)>0, 1, 'last');
end
CombinationsIndex_CombinationwithMaxNumberofSensors = max(numberofSensorsinCombinations);
SuccessfulSensorCombinaitons(:, CombinationsIndex_CombinationwithMaxNumberofSensors+1:end) = [];







% Save DPACoverageCombinationPriorityList
DPACoverageCombinationPriorityList_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPACoverageCombinationPriorityList.mat');


DPACoverageCombinationPriorityList_fullFileName = fullfile(ELS_dataSaveFolder, DPACoverageCombinationPriorityList_fileName);
save(DPACoverageCombinationPriorityList_fullFileName, 'DPACoverageCombinationPriorityList', 'SuccessfulSensorCombinaitons');
clear DPACoverageCombinationPriorityList_fileName DPACoverageCombinationPriorityList_fullFileName













xx = 1;       fprintf('\nESC Logic System: ESC sensor Coverage (DPACov) Combination Priority List :: Complete\n')    




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate ESC Logic System Information and Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ESCsensorInvolvement, DPACovCombination] = Function_ELSAnalysis_InfoStatGeneration_ST2(DPA_Zone, DPA_ID); %#ok<ASGLU>



















xx=1;       fprintf('\nESC Logic System: Information and Statistics Generation :: Complete\n')      % Completes gathering ESC sensor information

































end














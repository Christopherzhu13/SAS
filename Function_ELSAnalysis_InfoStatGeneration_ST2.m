%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_ELSAnalysis_InfoStatGeneration_ST2.m
%%%%
%%%%        - List of ESC Sensor Combination Indices
%%%%            - List of Indexes of combinations for which all the members sensors are Principal sensors
%%%%            - List of Indexes of combinations for which all the members sensors are Neighboring sensors (if any)
%%%%            - List of Indexes of combinations for which the members sensors consists of both Princiapal sensors AND Neighboring sensors
%%%%
%%%%        - List of Successfull ESC Sensor Combination (DPACov Combination) Indices
%%%%            - List of Indexes of DPACov Combinations for which all the members sensors are Principal sensors
%%%%            - List of Indexes of DPACov Combinations for which all the members sensors are Neighboring sensors (if any)
%%%%            - List of Indexes of DPACov Combinations for which the members sensors consists of both Princiapal sensors AND Neighboring sensors
%%%%
%%%%        - List of DPACov Combination with specific member Sensor
%%%%            - List of Indexes of combinations for which a specific ESC sensor is a member
%%%%            - List of Indexes of Principal combinations for which a specific Principal ESC sensor is member
%%%%
%%%%        - List of Available DPACov Combination with specific member Sensor unavailable
%%%%            - List of Indexes of combinations available for DPA Coverage if 'One' specific ESC sensor if 'OFF'
%%%%            - List of Indexes of Principal sensor combinations available for DPA Coverage if 'One' specific ESC sensor if 'OFF'
%%%%            - List of Indexes of combinations available for DPA Coverage if 'Two' specific ESC sensor if 'OFF'
%%%%            - List of Indexes of Principal sensor combinations available for DPA Coverage if 'Two' specific ESC sensor if 'OFF'
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
%%%%        - Generates and saves ESL Information and Statistics
%%%%            - ESCLogicSystem_InfoandStats.mat
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [ESCsensorInvolvement, DPACovCombination] = Function_ELSAnalysis_InfoStatGeneration_ST2(DPA_Zone, DPA_ID)



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
DPA_folderName = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));














xx=1;   % Completes Input initialization
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating filenames to load ESC Logic System data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELSFolder = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization';
addpath(ELSFolder)
addpath(genpath(ELSFolder))

ELS_dataSaveFolder =strcat(ELSFolder, '\Data_20190531\', DPA_folderName);

fileName_ESCsensorInfo = strcat(DPA_Name, '_ESCsensorInfo.mat');
fileName__DPACoverageCombinationPriorityList = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPACoverageCombinationPriorityList.mat');



fullFileName_ESCsensorInfo = fullfile(ELS_dataSaveFolder, fileName_ESCsensorInfo);
fullFileName_DPACoverageCombinationPriorityList = fullfile(ELS_dataSaveFolder, fileName__DPACoverageCombinationPriorityList);



load(fullFileName_ESCsensorInfo);                               % ESCsensorInfo: for each of the ESC sensors: Struct with Lon, Lat, Height, Azimuth, SiteID, SensorCategory, CoverageReliability
load(fullFileName_DPACoverageCombinationPriorityList);          % SuccessfulSensorCombinaitons and DPACoverageCombinationPriorityList: 
                                                                % DPACoverageCombinationPriorityList: Struct with Fields:
                                                                %       SensorCombination: cell{ for each successful combination, member sensor Indices according to their priority: [ESCsensorCombination]}
                                                                %       GPCoverageInfo: cell{ for each successful combination: [GPIndex, CombinationCoverageReliability]}



clear fileName_ESCsensorInfo fullFileName_ESCsensorInfo
clear fileName_DPACoverageCombinationPriorityList fullFileName_DPACoverageCombinationPriorityList















xx=1;   % Completes ESL data loading















%%% Analyze the ESC Sensor involvement in the DPA Coverage Performance of the ESC Logic System
% Find the number and indices of the ESC sensors involved in the ESC Logic System Analysis
escCount = size(ESCsensorInfo.SiteID,1);                                    % Total Number of ESC sensors considered in the ESL Analysis
AllSensorList = (1:escCount)';                                              % Indices of the ESC sensors used in the ESL Analysis

principalESCcount = sum(ESCsensorInfo.SensorCategory==1);                   % Number of Principal ESC sensors considered in the ESL Analysis
PrincipalSensorList = (1:principalESCcount)';                               % Indices of the Principal ESC sensors used in the ESL Analysis

NeighboringESCcount = sum(ESCsensorInfo.SensorCategory==0);                 % Number of Neighboring ESC sensors considered in the ESL Analysis
NeighboringSensorList = principalESCcount + (1:NeighboringESCcount)';       % Indices of the NEighboring ESC sensors used in the ESL Analysis




% Find the number and indices of the ESC sensors involved in the DPACov combinations
SuccessfulCombinaitonAllESCsensorList = unique(SuccessfulSensorCombinaitons(ismember(SuccessfulSensorCombinaitons, AllSensorList)));%#ok<NODEF>         % Indices of the ESC sensors involved in DPACov combinations

SuccessfulCombinaitonPrincipalESCsensorList = unique(SuccessfulSensorCombinaitons(ismember(SuccessfulSensorCombinaitons, PrincipalSensorList)));        % Indices of the Principal ESC sensors involved in DPACov combinations

SuccessfulCombinaitonNeighboringESCsensorList = unique(SuccessfulSensorCombinaitons(ismember(SuccessfulSensorCombinaitons, NeighboringSensorList)));    % Indices of the Neighboring ESC sensors involved in DPACov combinations



% Record the ESC sensor involvement information
ESCsensorInvolvement.AvailableSensorIndexList{1,1} = AllSensorList;
ESCsensorInvolvement.AvailableSensorIndexList{1,2} =  ESCsensorInfo.SiteID(AllSensorList);
ESCsensorInvolvement.AvailablePrincipalSensorIndexList{1,1} = PrincipalSensorList;
ESCsensorInvolvement.AvailablePrincipalSensorIndexList{1,2} =  ESCsensorInfo.SiteID(PrincipalSensorList);
ESCsensorInvolvement.AvailableNeighboringSensorIndexList{1,1} = NeighboringSensorList;
ESCsensorInvolvement.AvailableNeighboringSensorIndexList{1,2} =  ESCsensorInfo.SiteID(NeighboringSensorList);

ESCsensorInvolvement.SensorsinDPACovCombinations{1,1} = SuccessfulCombinaitonAllESCsensorList;
ESCsensorInvolvement.SensorsinDPACovCombinations{1,2} =  ESCsensorInfo.SiteID(SuccessfulCombinaitonAllESCsensorList);
ESCsensorInvolvement.PrincipalSensorsinDPACovCombinations{1,1} = SuccessfulCombinaitonPrincipalESCsensorList;
ESCsensorInvolvement.PrincipalSensorsinDPACovCombinations{1,2} =  ESCsensorInfo.SiteID(SuccessfulCombinaitonPrincipalESCsensorList);
ESCsensorInvolvement.NeighboringSensorsinDPACovCombinations{1,1} = SuccessfulCombinaitonNeighboringESCsensorList;
ESCsensorInvolvement.NeighboringSensorsinDPACovCombinations{1,2} =  ESCsensorInfo.SiteID(SuccessfulCombinaitonNeighboringESCsensorList);







% Save ESCsensorInvolvement
ESCsensorInvolvement_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_ESCsensorInvolvement.mat');
ESCsensorInvolvement_fullFileName = fullfile(ELS_dataSaveFolder, ESCsensorInvolvement_fileName);
save(ESCsensorInvolvement_fullFileName, 'ESCsensorInvolvement');
clear ESCsensorInvolvement_fileName ESCsensorInvolvement_fullFileName


















xx=1;       fprintf('\nESC Logic System: Info&Stat Collection on ESC sensor involvement in the DPA Coverage performance :: Complete\n')

















%%% Number of successful ESC sensor Combination with a) Principal Sensors only, b) Neighboring sensors only, and c) both Princial and Neighboring sensors
SuccessfulCombinationCount = size(DPACoverageCombinationPriorityList.SensorCombination,1);
AllsensorSuccessfulCombinationIndexList = (1:SuccessfulCombinationCount)';
PrincipalSensorSuccessfulCombinationIndexList = find(all(ismember(SuccessfulSensorCombinaitons, [PrincipalSensorList;0]),2));
NeighboringSensorSuccessfulCombinationIndexList = find(all(ismember(SuccessfulSensorCombinaitons, [NeighboringSensorList;0]),2));
CrossSensorSuccessfulCombinationIndexList = setdiff(AllsensorSuccessfulCombinationIndexList, ...
    [PrincipalSensorSuccessfulCombinationIndexList; NeighboringSensorSuccessfulCombinationIndexList]);


SuccessfulPrincipalSensorCombinaitons = SuccessfulSensorCombinaitons(PrincipalSensorSuccessfulCombinationIndexList, :); %#ok<IDISVAR>
SuccessfulNeighboringSensorCombinaitons = SuccessfulSensorCombinaitons(NeighboringSensorSuccessfulCombinationIndexList, :);
SuccessfulCrossSensorCombinaitons = SuccessfulSensorCombinaitons(CrossSensorSuccessfulCombinationIndexList, :);


DPACovCombination.SuccessfulAllSensorCombinations{1,1} = AllsensorSuccessfulCombinationIndexList;
DPACovCombination.SuccessfulAllSensorCombinations{1,2} = SuccessfulSensorCombinaitons;
DPACovCombination.SuccessfulPrincipalSensorCombinations{1,1} = PrincipalSensorSuccessfulCombinationIndexList;
DPACovCombination.SuccessfulPrincipalSensorCombinations{1,2} = SuccessfulPrincipalSensorCombinaitons;
DPACovCombination.SuccessfulNeighboringSensorCombinaitons{1,1} = NeighboringSensorSuccessfulCombinationIndexList;
DPACovCombination.SuccessfulNeighboringSensorCombinaitons{1,2} = SuccessfulNeighboringSensorCombinaitons;
DPACovCombination.SuccessfulCrossSensorCombinaitons{1,1} = CrossSensorSuccessfulCombinationIndexList;
DPACovCombination.SuccessfulCrossSensorCombinaitons{1,2} = SuccessfulCrossSensorCombinaitons;







% Save DPACovCombination
DPACovCombination_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPACovCombination.mat');
DPACovCombination_fullFileName = fullfile(ELS_dataSaveFolder, DPACovCombination_fileName);
save(DPACovCombination_fullFileName, 'DPACovCombination');
clear DPACovCombination_fileName DPACovCombination_fullFileName














xx=1;       fprintf('\nESC Logic System: Info&Stat Collection on DPA Coverage Combination classification based on sensor category :: Complete\n')















%%% Create the DPAov Combination table to be saved as .xlsx file
% Table for DPACov Combinations including all the ESC sensors involved
for combinationIndex = 1 : size(SuccessfulSensorCombinaitons,1)
    
    for sensorIndex = 1 : size(SuccessfulSensorCombinaitons,2)
        
        iterationSensorIndex = SuccessfulSensorCombinaitons(combinationIndex, sensorIndex);
        
        if iterationSensorIndex==0
            SuccessfulSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' ( - )'];
        else
            SuccessfulSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' (', ESCsensorInfo.SiteID{iterationSensorIndex}, ')'];
        end
        
    end
    
end


CombinationPriority = (1:size(SuccessfulSensorCombinaitons,1))';
ELSSensorCombinationSensorSiteIDTable = table(CombinationPriority, SuccessfulSensorCombinaitonsSiteIDs);
ELSSensorCombinationSensorSiteIDTable.Properties.VariableNames{1} = 'CombinationPriority';
ELSSensorCombinationSensorSiteIDTable.Properties.VariableNames{2} = 'ESCsensorSiteID';

% Save the ESC sensor combination Information (.xlsx file)
ELSSensorCombinationsSiteID_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_CombinationSensorSiteID.xlsx');
ELSSensorCombinationsSiteID_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorCombinationsSiteID_filename);

writetable(ELSSensorCombinationSensorSiteIDTable, ELSSensorCombinationsSiteID_fullFileName, 'Sheet', 1, 'Range','C5')





% NTIA software input generation: DPAname_ESC_Combo_Logic_Inputs.xlsx
for combinationIndex = 1 : size(SuccessfulSensorCombinaitons,1)    
    for sensorIndex = 1 : size(SuccessfulSensorCombinaitons,2)
        iterationSensorIndex = SuccessfulSensorCombinaitons(combinationIndex, sensorIndex);
        
        if iterationSensorIndex==0
            SuccessfulSensorCombinaitonsSiteIDs2{combinationIndex, sensorIndex} = ' ';
        else
            SuccessfulSensorCombinaitonsSiteIDs2{combinationIndex, sensorIndex} = ESCsensorInfo.SiteID{iterationSensorIndex};        
        end
    end    
end

if strcmp(DPA_Zone, 'east')
    DPAcovered_Zone = 'East';
elseif strcmp(DPA_Zone, 'west')
    DPAcovered_Zone = 'West';
end        
DPAcovered_Name = [DPAcovered_Zone, ' ', num2str(DPA_ID)];
DPACovered = repmat(DPAcovered_Name, size(SuccessfulSensorCombinaitons,1), 1);

if (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, [28, 29])) || (strcmp(DPA_Zone, 'west') && DPA_ID==16)
    ELSSensorCombinationSensorSiteIDTable2_cellArray{1,1} = DPACovered;
    ELSSensorCombinationSensorSiteIDTable2_cellArray{1,2} = SuccessfulSensorCombinaitonsSiteIDs2{1,1};
    ELSSensorCombinationSensorSiteIDTable2 = cell2table(ELSSensorCombinationSensorSiteIDTable2_cellArray);
else
    ELSSensorCombinationSensorSiteIDTable2 = table(DPACovered, SuccessfulSensorCombinaitonsSiteIDs2);
end
ELSSensorCombinationSensorSiteIDTable2.Properties.VariableNames{1} = 'DPACovered';
ELSSensorCombinationSensorSiteIDTable2.Properties.VariableNames{2} = 'SiteName';

% Save the ESC sensor combination Information (.xlsx file)
ELSSensorCombinationsSiteID_filename2 = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESC_Combo_Logic_Inputs.xlsx');
ELSSensorCombinationsSiteID_fullFileName2 = fullfile(ELS_dataSaveFolder, ELSSensorCombinationsSiteID_filename2);

writetable(ELSSensorCombinationSensorSiteIDTable2, ELSSensorCombinationsSiteID_fullFileName2, 'Sheet', 1, 'Range','A1')







% Table for Principal Sensor DPACov Combinations
if ~isempty(SuccessfulPrincipalSensorCombinaitons)
    for combinationIndex = 1 : size(SuccessfulPrincipalSensorCombinaitons,1)
        
        PrincipalSensorSuccessfulCombination_iteration = SuccessfulPrincipalSensorCombinaitons(combinationIndex, :);
        
        for sensorIndex = 1 : size(PrincipalSensorSuccessfulCombination_iteration,2)
            
            iterationSensorIndex = PrincipalSensorSuccessfulCombination_iteration(1, sensorIndex);
            
            if iterationSensorIndex==0
                SuccessfulPrincipalSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' ( - )'];
            else
                SuccessfulPrincipalSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' (', ESCsensorInfo.SiteID{iterationSensorIndex}, ')'];
            end
            
        end
        
    end
    
    
    PrincipalSensorCombinationPriority = PrincipalSensorSuccessfulCombinationIndexList;
    ELSPrincipalSensorCombinationSensorSiteIDTable = table(PrincipalSensorCombinationPriority, SuccessfulPrincipalSensorCombinaitonsSiteIDs);
    ELSPrincipalSensorCombinationSensorSiteIDTable.Properties.VariableNames{1} = 'PrincipalCombinationPriority';
    ELSPrincipalSensorCombinationSensorSiteIDTable.Properties.VariableNames{2} = 'PrincipalESCsensorSiteID';
    
    % Save the Principal ESC sensor combination Information (.xlsx file)
    ELSSensorCombinationsSiteID_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_CombinationSensorSiteID.xlsx');
    ELSSensorCombinationsSiteID_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorCombinationsSiteID_filename);
    
    writetable(ELSPrincipalSensorCombinationSensorSiteIDTable, ELSSensorCombinationsSiteID_fullFileName, 'Sheet', 2, 'Range','C5')
    
end







% Table for Neighboring Sensor DPACov Combinations
if ~isempty(SuccessfulNeighboringSensorCombinaitons)
    for combinationIndex = 1 : size(SuccessfulNeighboringSensorCombinaitons,1)
        
        NeighboringSensorSuccessfulCombination_iteration = SuccessfulNeighboringSensorCombinaitons(combinationIndex, :);
        
        for sensorIndex = 1 : size(NeighboringSensorSuccessfulCombination_iteration,2)
            
            iterationSensorIndex = NeighboringSensorSuccessfulCombination_iteration(1, sensorIndex);
            
            if iterationSensorIndex==0
                SuccessfulNeighboringSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' ( - )'];
            else
                SuccessfulNeighboringSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' (', ESCsensorInfo.SiteID{iterationSensorIndex}, ')'];
            end
            
        end
        
    end
    
    
    NeighboringSensorCombinationPriority = NeighboringSensorSuccessfulCombinationIndexList;
    ELSNeighboringSensorCombinationSensorSiteIDTable = table(NeighboringSensorCombinationPriority, SuccessfulNeighboringSensorCombinaitonsSiteIDs);
    ELSNeighboringSensorCombinationSensorSiteIDTable.Properties.VariableNames{1} = 'NeighboringCombinationPriority';
    ELSNeighboringSensorCombinationSensorSiteIDTable.Properties.VariableNames{2} = 'NeighboringESCsensorSiteID';
    
    % Save the Principal ESC sensor combination Information (.xlsx file)
    ELSSensorCombinationsSiteID_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_CombinationSensorSiteID.xlsx');
    ELSSensorCombinationsSiteID_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorCombinationsSiteID_filename);
    
    NeighboringCombination_RangeValue = strcat('C', num2str(5+size(SuccessfulPrincipalSensorCombinaitons,1)+10));    
    writetable(ELSNeighboringSensorCombinationSensorSiteIDTable, ELSSensorCombinationsSiteID_fullFileName, 'Sheet', 2, 'Range', NeighboringCombination_RangeValue)
    
end








% Table for Cross Sensor DPACov Combinations
if ~isempty(SuccessfulCrossSensorCombinaitons)
    for combinationIndex = 1 : size(SuccessfulCrossSensorCombinaitons,1)
        
        CrossSensorSuccessfulCombination_iteration = SuccessfulCrossSensorCombinaitons(combinationIndex, :);
        
        for sensorIndex = 1 : size(CrossSensorSuccessfulCombination_iteration,2)
            
            iterationSensorIndex = CrossSensorSuccessfulCombination_iteration(1, sensorIndex);
            
            if iterationSensorIndex==0
                SuccessfulCrossSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' ( - )'];
            else
                SuccessfulCrossSensorCombinaitonsSiteIDs{combinationIndex, sensorIndex} = [num2str(iterationSensorIndex), ' (', ESCsensorInfo.SiteID{iterationSensorIndex}, ')'];
            end
            
        end
        
    end
    
    
    CrossSensorCombinationPriority = CrossSensorSuccessfulCombinationIndexList;
    ELSCrossSensorCombinationSensorSiteIDTable = table(CrossSensorCombinationPriority, SuccessfulCrossSensorCombinaitonsSiteIDs);
    ELSCrossSensorCombinationSensorSiteIDTable.Properties.VariableNames{1} = 'CrossCombinationPriority';
    ELSCrossSensorCombinationSensorSiteIDTable.Properties.VariableNames{2} = 'CrossESCsensorSiteID';
    
    % Save the Principal ESC sensor combination Information (.xlsx file)
    ELSSensorCombinationsSiteID_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_CombinationSensorSiteID.xlsx');
    ELSSensorCombinationsSiteID_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorCombinationsSiteID_filename);
    
    CrossCombination_RangeValue = strcat('C', num2str(5+size(SuccessfulPrincipalSensorCombinaitons,1)+size(SuccessfulPrincipalSensorCombinaitons,1)+20));    
    writetable(ELSCrossSensorCombinationSensorSiteIDTable, ELSSensorCombinationsSiteID_fullFileName, 'Sheet', 2, 'Range', CrossCombination_RangeValue)
    
end





















xx=1;       fprintf('\nESC Logic System: Creating the DPAov Combination table and saving them as .xlsx file :: Complete\n')






















% Frequency and Significance of individual ESC sensors appearing in the DPA Coverage Combinations
DPACovCombination_PriorityFactor = (SuccessfulCombinationCount:-1:1)';

SuccessfulCombinationSensorFrequency = zeros(escCount, 1);
SuccessfulCombinationSensorSignificance = zeros(escCount, 1);
for successfulCombinationIndex = 1 : size(SuccessfulSensorCombinaitons,1)
    
    currentCombination = SuccessfulSensorCombinaitons(successfulCombinationIndex, :);
    currentCombination(currentCombination==0) = [];    
    
    SuccessfulCombinationSensorFrequency(currentCombination,1) = SuccessfulCombinationSensorFrequency(currentCombination,1) + 1;
    
    combinationPriorityIndex = DPACovCombination.SuccessfulAllSensorCombinations{1,1}(successfulCombinationIndex,1);
    SuccessfulCombinationSensorSignificance(currentCombination,1) = SuccessfulCombinationSensorSignificance(currentCombination,1) ...
        + DPACovCombination_PriorityFactor(combinationPriorityIndex,1);
    
end


SuccessfulCombinationSensorSignificance = (SuccessfulCombinationSensorSignificance ./ max(SuccessfulCombinationSensorSignificance)) * 10;










%%% Generate and save the table for DPA Coverage Combination Sensor Frequency to be saved as .xlsx file
[SuccessfulCombinationSensorFrequency_sorted, SuccessfulCombinationSensorFrequency_sortedIndex] = sort(SuccessfulCombinationSensorFrequency, 'descend');
SuccessfulCombinationSensorSiteID_sorted = ESCsensorInfo.SiteID(SuccessfulCombinationSensorFrequency_sortedIndex);

SuccessfulCombinationSensorSensorCategory_sortedTemp = ESCsensorInfo.SensorCategory(SuccessfulCombinationSensorFrequency_sortedIndex);
SuccessfulCombinationSensorSensorCategory_sorted = strings(escCount,1);
for escIndex = 1 : escCount
    if SuccessfulCombinationSensorSensorCategory_sortedTemp(escIndex,1) == 1
        SuccessfulCombinationSensorSensorCategory_sorted(escIndex,1) = 'Principal';
    else
        SuccessfulCombinationSensorSensorCategory_sorted(escIndex,1) = 'Neighboring';
    end
end


ELSSensorFrequencyTable = table(SuccessfulCombinationSensorSensorCategory_sorted, SuccessfulCombinationSensorSiteID_sorted, SuccessfulCombinationSensorFrequency_sorted); 
ELSSensorFrequencyTable.Properties.VariableNames{1} = 'SensorCategory';
ELSSensorFrequencyTable.Properties.VariableNames{2} = 'SensorSiteID';
ELSSensorFrequencyTable.Properties.VariableNames{3} = 'AssociatedCombinationNumber';

% Save the DPA Coverage Combination Sensor Frequency Information (.xlsx file)
ELSSensorFrequencyTable_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_SensorFrequency.xlsx');
ELSSensorFrequencyTable_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorFrequencyTable_filename);

writetable(ELSSensorFrequencyTable, ELSSensorFrequencyTable_fullFileName, 'Sheet',1, 'Range','C5')
clear ELSSensorFrequencyTable_filename ELSSensorFrequencyTable_fullFileName












%%% Generate and save the table for DPA Coverage Combination Sensor Significance to be saved as .xlsx file
[SuccessfulCombinationSensorSignificance_sorted, SuccessfulCombinationSensorSignificance_sortedIndex] = sort(SuccessfulCombinationSensorSignificance, 'descend');
SuccessfulCombinationSensorFrequency_sorted = SuccessfulCombinationSensorFrequency(SuccessfulCombinationSensorSignificance_sortedIndex,1);
SuccessfulCombinationSensorSiteID_sorted = ESCsensorInfo.SiteID(SuccessfulCombinationSensorSignificance_sortedIndex);

SuccessfulCombinationSensorSensorCategory_sortedTemp = ESCsensorInfo.SensorCategory(SuccessfulCombinationSensorSignificance_sortedIndex);
SuccessfulCombinationSensorSensorCategory_sorted = strings(escCount,1);
for escIndex = 1 : escCount
    if SuccessfulCombinationSensorSensorCategory_sortedTemp(escIndex,1) == 1
        SuccessfulCombinationSensorSensorCategory_sorted(escIndex,1) = 'Principal';
    else
        SuccessfulCombinationSensorSensorCategory_sorted(escIndex,1) = 'Neighboring';
    end
end

format short
ELSSensorSignificanceTable = table(SuccessfulCombinationSensorSensorCategory_sorted, ...
    SuccessfulCombinationSensorSiteID_sorted, SuccessfulCombinationSensorFrequency_sorted, floor(SuccessfulCombinationSensorSignificance_sorted)); 
ELSSensorSignificanceTable.Properties.VariableNames{1} = 'SensorCategory';
ELSSensorSignificanceTable.Properties.VariableNames{2} = 'SensorSiteID';
ELSSensorSignificanceTable.Properties.VariableNames{3} = 'AssociatedCombinationNumber';
ELSSensorSignificanceTable.Properties.VariableNames{4} = 'SensorSignificance_max10';


% Save the DPA Coverage Combination Sensor Frequency Information (.xlsx file)
ELSSensorSignificanceTable_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_SensorFrequency.xlsx');
ELSSensorSignificanceTable_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorSignificanceTable_filename);

writetable(ELSSensorSignificanceTable, ELSSensorSignificanceTable_fullFileName, 'Sheet',2, 'Range','C5')
clear ELSSensorSignificanceTablee_filename ELSSensorSignificanceTable_fullFileName












%%% Generate and save the table for DPA Coverage Combination Sensor Significance to be saved as .xlsx file
ELSSensorSignificanceTable_Sheet3 = sortrows(ELSSensorSignificanceTable, 1, 'descend');

% Save the DPA Coverage Combination Sensor Frequency Information (.xlsx file)
ELSSensorSignificanceTable_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCLogicSystem_SensorFrequency.xlsx');
ELSSensorSignificanceTable_fullFileName = fullfile(ELS_dataSaveFolder, ELSSensorSignificanceTable_filename);

writetable(ELSSensorSignificanceTable_Sheet3, ELSSensorSignificanceTable_fullFileName, 'Sheet', 3, 'Range', 'C5')
clear ELSSensorSignificanceTablee_filename ELSSensorSignificanceTable_fullFileName





















xx=1;       fprintf('\nESC Logic System: Generate and save the table for DPA Coverage Combination Sensor Frequency to be saved as .xlsx file :: Complete\n')























% Frequency of Principal ESC sensors appearing in the DPA Coverage  Principal Combinations
% SuccessfulPrincipalCombinationSensorFrequency = zeros(principalESCcount, 1);
% SuccessfulPrincipalCombinationSensorSignificance = zeros(escCount, 1);
% for successfulPrincipalCombinationIndex = 1 : size(SuccessfulPrincipalSensorCombinaitons,1)
%     
%     currentPrincipalCombination = SuccessfulPrincipalSensorCombinaitons(successfulPrincipalCombinationIndex, :);
%     currentPrincipalCombination(currentPrincipalCombination==0) = [];    
%     
%     SuccessfulPrincipalCombinationSensorFrequency(currentPrincipalCombination,1) = SuccessfulPrincipalCombinationSensorFrequency(currentPrincipalCombination,1) + 1;
%     
%     combinationPriorityIndex = DPACovCombination.SuccessfulPrincipalSensorCombinations{1,1}(successfulPrincipalCombinationIndex,1);
%     SuccessfulPrincipalCombinationSensorSignificance(currentPrincipalCombination,1) = SuccessfulPrincipalCombinationSensorSignificance(currentPrincipalCombination,1) ...
%         + DPACovCombination_PriorityFactor(combinationPriorityIndex,1);
%     
% end






















xx=1;



































end















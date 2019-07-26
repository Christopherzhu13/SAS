%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% ESLAnalysis_CombinationCoveragePerformance.m
%%%% ESC Sensor Logic Analysis Platform
%%%%    1. Provides sensor combination coverage performance
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear all; clc; %#ok<CLALL>
format long;



% east_dpa_1.geojson (1 to 23); 
% east_dpa_bath.geojson; 
% east_dpa_mayport.geojson; 
% east_dpa_norfolk.geojson; 
% east_dpa_pascagoula.geojson;
% west_dpa_1.geojson (1 to 14); 
% west_dpa_bremerton.geojson
% alaska_dpa_1.geojson (1 to 38)
% hawaii_dpa_1.geojson (1 to 9)
% puerto_rico_ver1_dpa_1.geojson (1 to 5); 
% puerto_rico_ver2_dpa_1.geojson (1 to 5); 
% puerto_rico_DPAs-9-29-2017-Version-1.geojson; 
% puerto_rico_DPAs-9-29-2017-Version-2.geojson

% NOT DONE:     east_dpa_bath.geojson;                  No ESC Info









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taking user inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Please select the DPA Zone and ID for the Coverage Analysis. \n\n')
fprintf('\nOptions for DPA Zone: east, west \n\n')

prompt = 'Please enter the DPA zone: ';
DPA_Zone = input(prompt, 's');
if isempty(DPA_Zone)
    DPA_Zone = 'west';
end








fprintf('\n\nOptions for DPA ID: \n')
fprintf('For DPA Zone = east, DPA ID options are: 1 - 23, 24 for Mayport, 25 for Norfolk, 26 for Pascagoula, 27 for Beth \n')
fprintf('For DPA Zone = west, DPA ID options are: 1 - 14, 15 for bremerton \n')
fprintf('Not Functional: \n')
fprintf('   East DPA Beth:  no ESC information available \n')


prompt = '\nPlease enter the DPA ID: ';
DPA_ID = input(prompt);
if isempty(DPA_ID)
    DPA_ID = 3;
end



















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
DPAname = strcat(DPA_Zone, '_', num2str(DPA_IDstr));
DPAfolderName = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
ESCInfoDPAname = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));











xx=1;   % Completes Input initialization










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating filenames to load DPA Coverage Analysis data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELSFolder = 'C:\Users\czhu\Dropbox (Federated Wireless)\Share with Chenhan\Azimuth Optimization';
addpath(ELSFolder)
addpath(genpath(ELSFolder))

ELS_dataFolder =strcat(ELSFolder, '\Data_20190625\', DPAfolderName);

fileName_DPABoundary = strcat(DPAname, '_DPA_Boundary.mat');
fileName_DPAinnerBoundary = strcat(DPAname, '_DPA_innerBoundary.mat');
fileName_DPAgridInfo = strcat(DPAname, '_DPA_gridInfo.mat');
fileName_ESCInfo = strcat(DPAname, '_ESCsensorInfo.mat');
fileName_GridPointDetectionandCoverageInfo = strcat(DPAname, '_GridPointDetectionandCoverageInfo.mat');
fileName_DPACoverageCombinationPriorityList = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPACoverageCombinationPriorityList.mat');
fileName_ESCsensorInvolvement = strcat(DPAname, '_ESCsensorInvolvement.mat');
% fileName_DPACovCombination = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_DPACovCombination.mat');


fullFileName_DPABoundary = fullfile(ELS_dataFolder, fileName_DPABoundary);
fullFileName_DPAinnerBoundary = fullfile(ELS_dataFolder, fileName_DPAinnerBoundary);
fullFileName_DPAgridInfo = fullfile(ELS_dataFolder, fileName_DPAgridInfo);
fullFileName_ESCInfo = fullfile(ELS_dataFolder, fileName_ESCInfo);
fullFileName_GridPointDetectionandCoverageInfo = fullfile(ELS_dataFolder, fileName_GridPointDetectionandCoverageInfo);
fullFileName_DPACoverageCombinationPriorityList = fullfile(ELS_dataFolder, fileName_DPACoverageCombinationPriorityList);
fullFileName_ESCsensorInvolvement = fullfile(ELS_dataFolder, fileName_ESCsensorInvolvement);
% fullFileName_DPACovCombination = fullfile(ELS_dataFolder, fileName_DPACovCombination);



load(fullFileName_DPABoundary);                             % DPA_Boundary: for each of the Boundary Points: [Lon, Lat, EdgeIndex, VertexIndex]
load(fullFileName_DPAinnerBoundary);                        % DPA_InnerBoundary: for each of the Inner Boundary Points: [Lon, Lat, InnerBoundaryPartIndex]
load(fullFileName_DPAgridInfo);                             % DPA_gridInfo: for each of the DPA Grid Points: [Lon, Lat, GPIndex, inside75kmPolygonGPIndex]
load(fullFileName_ESCInfo);                                 % ESCsensorInfo: for each of the ESC sensors: Struct with Lon, Lat, Height, Azimuth, SiteID, SensorCategory, CoverageReliability
load(fullFileName_GridPointDetectionandCoverageInfo);       % GridPointDetectionInfo: for each of the ESC sensors: Cell{ for each of the GP under detection: [GPIndex, DetectionReliability]}
load(fullFileName_DPACoverageCombinationPriorityList);          % DPACoverageCombinationPriorityList: struct.SensorCombination, struct.GPCoverageInfo
                                                                % struct.SensorCombination: cell{ for each successful combination, member sensor Indices according to their priority: [ESCsensorCombination]}
                                                                % struct.GPCoverageInfo: cell{ for each successful combination: [GPIndex, CombinationCoverageReliability]}
load(fullFileName_ESCsensorInvolvement)
% load(fullFileName_DPACovCombination)


clear fileName_DPABoundary fullFileName_DPABoundary
clear fileName_DPAinnerBoundary fullFileName_DPAinnerBoundary
clear fileName_DPAgridInfo fullFileName_DPAgridInfo
clear fileName_ESCInfo fullFileName_ESCInfo
clear fileName_GridPointDetectionandCoverageInfo fullFileName_GridPointDetectionandCoverageInfo
clear fileName_DPACoverageCombinationPriorityList fullFileName_DPACoverageCombinationPriorityList
clear fileName_ESCsensorInvolvement fullFileName_ESCsensorInvolvement
% clear fileName_DPACovCombination fullFileName_DPACovCombination




%%%% Pre-procesing the loaded ESC Sensor Logic data
%   For each boundary point: DPA_Boundary = [Lon, Lat, EdgeIndex, VertexIndex]

% Get DPA Boundary and Island Boundary information
if strcmp(DPA_Zone, 'west') && DPA_ID == 13
      
    DPA_Boundary = DPA_Boundary_cell{1};
    DPA_Boundary_Island_1 = DPA_Boundary_cell{2};
    DPA_Boundary_Island_2 = DPA_Boundary_cell{3};
    DPA_Boundary_Island_3 = DPA_Boundary_cell{4};
    
end








% Get the useful DPA Coverge Analysis varibales
escCount = length(ESCsensorInfo.SiteID);         escTally = ones(escCount, 1);
SelectedCombinationMax = 20;
if size(DPACoverageCombinationPriorityList.SensorCombination,1)<SelectedCombinationMax
    SelectedCombinationMax = size(DPACoverageCombinationPriorityList.SensorCombination,1);
end
if strcmp(DPA_Zone, 'east') && DPA_ID==8
    escSensorsInvolvedIndex = [];
    for selectedCombinationIndex = 1 : SelectedCombinationMax
        escSensorsInvolvedIndex = [escSensorsInvolvedIndex; DPACoverageCombinationPriorityList.SensorCombination{selectedCombinationIndex}];
    end
    
    escSensorsInvolvedIndex = unique(escSensorsInvolvedIndex);
    escCountInvolved = size(escSensorsInvolvedIndex,1);
    
else
    escCountInvolved = size(ESCsensorInvolvement.SensorsinDPACovCombinations{1},1);
end    

escTallyInvolved = ones(escCountInvolved, 1);
DPA_gridInfo_Temp = DPA_gridInfo(1:1:end,:);
gridPointCount = size(DPA_gridInfo_Temp, 1);     gridPointTally = ones(gridPointCount, 1);
DPA_gridLon = DPA_gridInfo_Temp(:, 1);           DPA_gridLat = DPA_gridInfo_Temp(:, 2);







xx=1;










% DPA_Boundary
% DPA_innerBoundary
% DPA_gridInfo
% ESCsensorInfo
% GridPointDetectionInfo
% GridPointCoverageInfo
% DPACoverageCombinationPriorityList
% ESCsensorInvolvement
% DPACovCombination








fprintf('\n\nOptions for ESC Sensor Combination: \n')
fprintf(strcat('Number of sensor combination satisfying DPA Coverage Requirements: ', num2str(size(DPACoverageCombinationPriorityList.SensorCombination,1)), '\n'))

prompt = '\nPlease enter the Sensor combination Index: ';
SensorCombinationIndex = input(prompt);
if isempty(SensorCombinationIndex)
    SensorCombinationIndex = 1;
end





















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coverage Info and Result Text Box Positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coverage Info Text Box Position
TxtBox1_E_TL = [1, 2, 3, 4, 9, 18, 19, 24, 25, 26];                 TxtBox1_W_TL = [1, 10, 11];         
TxtBox1_E_BL = [13, 14, 15, 17, 22];                                        TxtBox1_W_BL = [8, 9, 12, 13, 15];
TxtBox1_E_BR = [6, 16];                                                 TxtBox1_W_BR = [5, 14];
TxtBox1_E_TR = [5, 7, 8, 10, 11, 12, 20, 21, 23];                       TxtBox1_W_TR = [2, 3, 4, 6, 7];

if (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox1_E_TL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox1_W_TL))
    TB1_Dimension = [0.13 0.88 0.195 0.045];        TB1_HAlignment = 'left';
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox1_E_BL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox1_W_BL))
    TB1_Dimension = [0.13 0.11 0.195 0.045];        TB1_HAlignment = 'left';                
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox1_E_BR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox1_W_BR))
    TB1_Dimension = [0.71 0.11 0.195 0.045];        TB1_HAlignment = 'right';
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox1_E_TR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox1_W_TR))
    TB1_Dimension = [0.71 0.88 0.195 0.045];        TB1_HAlignment = 'right';
end     % If ends making sure that the TextBoxes are correctly positioned




% With ECSInfo_NperLine ESC information per line, setting up all the Text Box positions
ECSInfo_NperLine = 10;
% strN_nrows = ceil(escCountInvolved/ECSInfo_NperLine);
strN_nrows = 1;
if TB1_Dimension(2) > 0.35      % Make sure the Coverage information TextBox is positioned correctly    
    TB1_Dimension_All = TB1_Dimension - [0 0.03 0 0] .* (0:strN_nrows+1)';        
elseif TB1_Dimension(2) < 0.35        
    TB1_Dimension_All = TB1_Dimension + [0 0.03 0 0] .* (strN_nrows+1:-1:0)';        
end     % If ends making sure that the TextBoxes are positioned correctly










% Coverage Result Text Box Position
TxtBox2_E_TL = [5, 6, 7, 8, 10, 11, 12];                                TxtBox2_W_TL = [2, 3, 4, 5, 6, 7, 8];         
TxtBox2_E_BL = [20, 21, 23, 24];                                        TxtBox2_W_BL = 15;
TxtBox2_E_BR = [2, 3, 4, 13, 14, 15, 22, 25, 26];                       TxtBox2_W_BR = 0;
TxtBox2_E_TR = [1, 9, 16, 17, 18, 19];                                  TxtBox2_W_TR = [1, 9, 10, 11, 12, 13, 14];

if (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox2_E_TL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox2_W_TL))
    TB2_Dimension = [0.13 0.88 0.195 0.045];        TB2_HAlignment = 'left';
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox2_E_BL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox2_W_BL))
    TB2_Dimension = [0.13 0.11 0.195 0.045];        TB2_HAlignment = 'left';                
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox2_E_BR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox2_W_BR))
    TB2_Dimension = [0.71 0.11 0.195 0.045];        TB2_HAlignment = 'right';
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox2_E_TR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox2_W_TR))
    TB2_Dimension = [0.71 0.88 0.195 0.045];        TB2_HAlignment = 'right';
end     % If ends making sure that the TextBoxes are correctly positioned



if TB2_Dimension(2) > 0.35      % Make sure the Coverage information TextBox is positioned correctly    
    TB2_Dimension_All = TB2_Dimension - [0 0.03 0 0] .* (0:2)';        
elseif TB2_Dimension(2) < 0.35        
    TB2_Dimension_All = TB2_Dimension + [0 0.03 0 0] .* (2:-1:0)';        
end     % If ends making sure that the TextBoxes are positioned correctly











% DPA Coverage ESC sensor combination Text Box Position
% TxtBox3_E_TL = [15, 22];                                                            TxtBox3_W_TL = [9, 12, 13];         
% TxtBox3_E_BL = [14, 17, 19, 26];                                                    TxtBox3_W_BL = [1, 2, 3, 5, 6, 7, 10, 11, 14, 15];
% TxtBox3_E_BR = [1, 5, 8, 7, 9, 10, 11, 12, 13, 16, 18, 20, 21, 23, 24];             TxtBox3_W_BR = 4;
% TxtBox3_E_TR = [2, 3, 4, 6, 25];                                                    TxtBox3_W_TR = 8;
% 
% if (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox3_E_TL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox3_W_TL))
%     TB3_Dimension = [0.13 0.88 0.195 0.045];           TB3_HAlignment = 'left';         TB3_PositionInfo = [1, 2, 1, -1];
% elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox3_E_BL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox3_W_BL))
%     TB3_Dimension = [0.13 0.11 0.195 0.045];           TB3_HAlignment = 'left';         TB3_PositionInfo = [1, 1, 1, 1];
% elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox3_E_BR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox3_W_BR))
%     TB3_Dimension = [0.71 0.11 0.195 0.045];           TB3_HAlignment = 'right';         TB3_PositionInfo = [2, 1, -1, 1];
% elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, TxtBox3_E_TR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, TxtBox3_W_TR))
%     TB3_Dimension = [0.71 0.88 0.195 0.045];           TB3_HAlignment = 'right';         TB3_PositionInfo = [2, 2, -1, -1];
% end     % If ends making sure that the TextBoxes are correctly positioned
% 
% 
% 
% if TB3_Dimension(2) > 0.35      % Make sure the Coverage information TextBox is positioned correctly    
%     TB3_Dimension_All = TB3_Dimension - [0 0.047 0 0] .* (0:5)';
% elseif TB3_Dimension(2) < 0.35
%     TB3_Dimension_All = TB3_Dimension + [0 0.047 0 0] .* (5:-1:0)';
% end     % If ends making sure that the TextBoxes are positioned correctly












% Mode Selection Box Position
mode_E_TL = [13, 14, 20, 21, 23];                                           mode_W_TL = 14;         
mode_E_BL = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 18, 25];            mode_W_BL = 4;
mode_E_BR = [17, 19];                                                       mode_W_BR = [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13];
mode_E_TR = [15, 22, 24, 26];                                               mode_W_TR = [5, 15];


if (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, mode_E_TL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, mode_W_TL))
    mode_dim1 = [0.13 0.86 0.11 0.045];            MS_HAlignment = 'left';             MS1_PositionInfo = [1, 2, 1, -1];
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, mode_E_BL)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, mode_W_BL))
    mode_dim1 = [0.13 0.11 0.11 0.045];            MS_HAlignment = 'left';             MS1_PositionInfo = [1, 1, 1, 1];
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, mode_E_BR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, mode_W_BR))
    mode_dim1 = [0.795 0.11 0.11 0.045];            MS_HAlignment = 'right';            MS1_PositionInfo = [2, 1, -1, 1];
elseif (strcmp(DPA_Zone, 'east') && ismember(DPA_ID, mode_E_TR)) || (strcmp(DPA_Zone, 'west') && ismember(DPA_ID, mode_W_TR))
    mode_dim1 = [0.795 0.86 0.11 0.045];            MS_HAlignment = 'right';            MS1_PositionInfo = [2, 2, -1, -1];
end     % If ends making sure that the Mode Selection Boxes are correctly positioned



% Setting up all the 3 Legend Box positions
if mode_dim1(2) > 0.5      % Make sure the Coverage information TextBox is positioned correctly    
    mode_dim1_all = mode_dim1 - [0 0.03 0 0] .* (0:2)';        
elseif mode_dim1(2) < 0.5        
    mode_dim1_all = mode_dim1 + [0 0.03 0 0] .* (2:-1:0)';        
end  

% mode_dim1_all(1,2) = mode_dim1_all(1,2) + 0.02;




























































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Presenting the result:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CoverageReliabilityThreshold95 = 0.95;           % Reliability threshold of 0.95 for DPA Coverage percentage calculation
CoverageReliabilityThreshold50 = 0;            % Reliability threshold of 0.50 for DPA Coverage percentage calculation
gridPointTally95 = gridPointTally;
gridPointTally50 = gridPointTally;



target_gridIndex_init = 500;

DPA_Boundary_color = [0.64 0.08 0.18];
DPA_Boundary_Islandcolor = [0.64 0.08 0.18];
DPA_Boundary_innerBoundarycolor = [0.64 0.08 0.18];

Rge9_ESCcolor = [0 0.5 0];
Rlt9_ESCcolor = [0, 0, 1];
Rbase_ESCcolor = [0 0 0];
Principal_ESCcolor = [0.3 0.3 0.3];
Neighboring_ESCcolor = [0.5 0.5 0.5];
Highlight_ESCcolor = [0.64 0.08 0.18];
Deactivated_ESCcolor = [1 0 0];

Rge95_GPcolor = [0 1 0];
Rge5lt95_GPCcolor = [0 0 1];
Rlt5_GPcolor = [1 0 0];
Rbase_GPcolor = [0 0 0];
Default_GPcolor = [0.6, 0.6, 0.6];


Statement_TB1color = [0.64 0.08 0.18];
Info_TB1color = [0 0 0.7];


Statement_TB2color = [0.64 0.08 0.18];
Info_TB2color = [0 0 0.7];
pass_TB2color = [0 0.5 0];
fail_TB2color = [1 0 0];


Statement_TB3color = [0.64 0.08 0.18];
Info_TB3color = [0 0 0];
HPSC_TB3color = [0 0 1];


Choice1_MScolor = [0 0 0.7];
Choice2_MScolor = [0 0 1];




GPMarkerSize = 0.5;








% Setting up ESC Location Text offset
if strcmp(DPA_Zone, 'east') && strcmp(DPA_ID, 'Norfolk')
    escLocTxt_offsetX = 0.03;
    escLocTxt_offsetY = 0.01;
    
elseif strcmp(DPA_Zone, 'west') && DPA_ID == 12
    escLocTxt_offsetX = 0.30;
    escLocTxt_offsetY = 0.05;
else
    escLocTxt_offsetX = 0.05;
    escLocTxt_offsetY = 0.05;
end




















%%% Start the DPA Protection Analysis visualization
hf = figure;        % maximize(hf);
hold on












%%% Plot the Boundary points for the selected DPA
DPABoundary_ploth = plot(DPA_Boundary(:,1), DPA_Boundary(:,2), '--', 'color', DPA_Boundary_color, 'linewidth', 1);


% For DPA west-12 and west-13 plot the Island Boundaries
if strcmp(DPA_Zone, 'west') && DPA_ID == 13
    
    DPAIslandBoundary_ploth(1) = plot(DPA_Boundary_Island_1(:,1), DPA_Boundary_Island_1(:,2), '--', 'color', DPA_Boundary_Islandcolor, 'linewidth', 1);
    DPAIslandBoundary_ploth(2) = plot(DPA_Boundary_Island_2(:,1), DPA_Boundary_Island_2(:,2), '--', 'color', DPA_Boundary_Islandcolor, 'linewidth', 1);
    DPAIslandBoundary_ploth(3) = plot(DPA_Boundary_Island_3(:,1), DPA_Boundary_Island_3(:,2), '--', 'color', DPA_Boundary_Islandcolor, 'linewidth', 1);
    
end














%%% Plot the Inner Boundary for the slected DPA
if (strcmp(DPA_Zone, 'east') && le(DPA_ID, 26)) || (strcmp(DPA_Zone, 'west') && le(DPA_ID, 14))
    
    
    inner_Boundary_ploth = gobjects(1, size(unique(DPA_innerBoundary(:,3)),1));
    if any(DPA_innerBoundary(:,3))
        
        innerBoundaryPart = unique(DPA_innerBoundary(:,3));
        
        DPA_innerBoundaryPartstartIndex = 1;
        for innerBoundaryPartIndex = 1 : size(innerBoundaryPart,1)
            
            DPA_innerBoundaryPartendIndex = find(DPA_innerBoundary(:,3) == innerBoundaryPartIndex, 1, 'last');
            
            inner_Boundary_ploth(innerBoundaryPartIndex) = plot(DPA_innerBoundary(DPA_innerBoundaryPartstartIndex:DPA_innerBoundaryPartendIndex,1), ...
                DPA_innerBoundary(DPA_innerBoundaryPartstartIndex:DPA_innerBoundaryPartendIndex,2), '--', 'color', DPA_Boundary_innerBoundarycolor, 'linewidth', 1);
            
            DPA_innerBoundaryPartstartIndex = DPA_innerBoundaryPartendIndex + 1;
        end
        
    else
        inner_Boundary_ploth = plot(DPA_innerBoundary(:,1), DPA_innerBoundary(:,2), '--', 'color', DPA_Boundary_innerBoundarycolor, 'linewidth', 1);
    end
    
end




















%%% Print the ESC sensors associated with the selected DPA
SelectedESCsensorCombination = DPACoverageCombinationPriorityList.SensorCombination{SensorCombinationIndex};

escLoc_ploth = gobjects(1, escCountInvolved);
escLocTxt_ploth = gobjects(1, escCountInvolved);
for escIndexInvolved = 1 : escCountInvolved         % For each of the ESCs plot the initial ESC information
    
    
    if strcmp(DPA_Zone, 'east') && DPA_ID==8
        escIndex = escSensorsInvolvedIndex(escIndexInvolved,1);
    else
        escIndex = ESCsensorInvolvement.SensorsinDPACovCombinations{1,1}(escIndexInvolved,1);
    end
    
    if ismember(escIndex, SelectedESCsensorCombination)
        ESCcolor = Highlight_ESCcolor;
        
    elseif ESCsensorInfo.SensorCategory(escIndex)==1
        ESCcolor = Principal_ESCcolor;
    else
        ESCcolor = Neighboring_ESCcolor;
    end
    escMarkerSize = 10;
    escLineWidth = 1.5;
    
    if strcmp(ESCsensorInfo.SiteID{escIndex}, 'NJ005')
        escLocTxt_offsetX_temp = escLocTxt_offsetX;
        escLocTxt_offsetY_temp = escLocTxt_offsetY;
        
        escLocTxt_offsetX = 0.05;
        escLocTxt_offsetY = -0.2;
        
    elseif strcmp(ESCsensorInfo.SiteID{escIndex}, 'RI003')
        escLocTxt_offsetX_temp = escLocTxt_offsetX;
        escLocTxt_offsetY_temp = escLocTxt_offsetY;
        
        escLocTxt_offsetX = -0.1;
        escLocTxt_offsetY = -0.1;
        
    elseif strcmp(ESCsensorInfo.SiteID{escIndex}, 'NY005')
        escLocTxt_offsetX_temp = escLocTxt_offsetX;
        escLocTxt_offsetY_temp = escLocTxt_offsetY;
        
        escLocTxt_offsetX = -0.1;
        escLocTxt_offsetY = -0.05;
        
    elseif strcmp(ESCsensorInfo.SiteID{escIndex}, 'NC006')
        escLocTxt_offsetX_temp = escLocTxt_offsetX;
        escLocTxt_offsetY_temp = escLocTxt_offsetY;
        
        escLocTxt_offsetX = 0.4;
        escLocTxt_offsetY = 0.05;
        
    elseif strcmp(ESCsensorInfo.SiteID{escIndex}, 'NC007') || strcmp(ESCsensorInfo.SiteID{escIndex}, 'NC008')
        escLocTxt_offsetX_temp = escLocTxt_offsetX;
        escLocTxt_offsetY_temp = escLocTxt_offsetY;
        
        escLocTxt_offsetX = 0.0;
        escLocTxt_offsetY = 0.05;
        
    end
    
    
        
    
    escLoc_ploth(escIndexInvolved) = plot(ESCsensorInfo.Longitude(escIndex), ESCsensorInfo.Latitude(escIndex), '*', ...
        'color', ESCcolor, 'markersize', escMarkerSize, 'linewidth', escLineWidth);
    
    txt = ['  ESC: ', ESCsensorInfo.SiteID{escIndex}];
    escLocTxt_ploth(escIndexInvolved) = text(ESCsensorInfo.Longitude(escIndex)+escLocTxt_offsetX, ESCsensorInfo.Latitude(escIndex)+escLocTxt_offsetY, txt, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'color', ESCcolor, 'FontWeight', 'bold');
    
    
    if strcmp(ESCsensorInfo.SiteID{escIndex}, 'NJ005') || strcmp(ESCsensorInfo.SiteID{escIndex}, 'RI003') || ...
            strcmp(ESCsensorInfo.SiteID{escIndex}, 'NY005') || strcmp(ESCsensorInfo.SiteID{escIndex}, 'NC007') || strcmp(ESCsensorInfo.SiteID{escIndex}, 'NC008')
        escLocTxt_offsetX = escLocTxt_offsetX_temp;
        escLocTxt_offsetY = escLocTxt_offsetY_temp;
    end
    
    
end     % For loop ends after printing all the ESC information



if strcmp(DPA_Zone, 'east') && DPA_ID==1
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1.1 xlimtemp(2)+0.8]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)]+0.1);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==2
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.6 xlimtemp(2)+0.1]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)+0.1 ylimtemp(2)+0.1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==3
    xlimtemp = xlim;
    xlim([xlimtemp(1)-2 xlimtemp(2)+2]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)+0.2]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==4
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.1 xlimtemp(2)+0.5]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==5
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1.1 xlimtemp(2)+1.2]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.2 ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==6
    xlimtemp = xlim;
    xlim([xlimtemp(1) xlimtemp(2)+0.2]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.3]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==7
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.4 xlimtemp(2)+1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==8
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.4 xlimtemp(2)+0.4]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.4]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==9
    xlimtemp = xlim;
    xlim([xlimtemp(1) xlimtemp(2)+0.4]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==10
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.1 xlimtemp(2)+0.2]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)+0.2 ylimtemp(2)-0.2]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==11
    xlimtemp = xlim;
    xlim([xlimtemp(1) xlimtemp(2)+0.8]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==12
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1.4 xlimtemp(2)+0.7]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==13
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.9 xlimtemp(2)+0.7]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==14
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.8 xlimtemp(2)+0.1]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==15
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.5 xlimtemp(2)+0.4]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==16
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.4 xlimtemp(2)]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.2 ylimtemp(2)+0.1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==17
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.4 xlimtemp(2)+1.2]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==18
    xlimtemp = xlim;
    xlim([xlimtemp(1) xlimtemp(2)+0.2]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==19
    xlimtemp = xlim;
    xlim([xlimtemp(1)+0.1 xlimtemp(2)+0.8]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)+0.1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==20
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1.5 xlimtemp(2)+1.5]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1)-0.1 ylimtemp(2)+0.1]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==21
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1 xlimtemp(2)+1.7]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.25]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==22
    xlimtemp = xlim;
    xlim([xlimtemp(1)-2 xlimtemp(2)+1.5]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.3]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==23
    xlimtemp = xlim;
    xlim([xlimtemp(1)+0.2 xlimtemp(2)-0.4]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.2]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==24
    xlimtemp = xlim;
    xlim([xlimtemp(1) xlimtemp(2)+1]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'east') && DPA_ID==26
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.8 xlimtemp(2)+1.3]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'west') && DPA_ID==5
    xlimtemp = xlim;
    xlim([xlimtemp(1)-0.8 xlimtemp(2)+1.3]);
    
    grid off;   grid on;    grid minor;
    
elseif strcmp(DPA_Zone, 'west') && DPA_ID==13
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1.3 xlimtemp(2)+1]);
    
    ylimtemp = ylim;
    ylim([ylimtemp(1) ylimtemp(2)+0.2]);
    
    grid off;   grid on;    grid minor;
    
else
    xlimtemp = xlim;
    xlim([xlimtemp(1)-1 xlimtemp(2)+1.0]);
    
    grid off;   grid on;    grid minor;
    
end


























%%% Color the DPA grid points according to the achieved coverage reliability
DPAgridPoint_ploth = gobjects(gridPointCount, 1);
GPTargetMarker_ploth = gobjects(gridPointCount, 1);     GPmarkerFlag = 0;
for DPA_gridIndex = 1 : gridPointCount          % Color each of the Grid Points according to the combined Grid Point RC
    
    delete(DPAgridPoint_ploth(DPA_gridIndex));
    DPAgridPoint_ploth(DPA_gridIndex) = plot(DPA_gridLon(DPA_gridIndex), DPA_gridLat(DPA_gridIndex), 'o', 'color', Default_GPcolor, ...
        'MarkerFaceColor', Default_GPcolor, 'MarkerSize', GPMarkerSize);
    
end



HPSC_GridPointReliability = DPACoverageCombinationPriorityList.GPcoverageInfo{SensorCombinationIndex}(:,2);
gridPointTally95_temp = gridPointTally95;
gridPointTally50_temp = gridPointTally50;

for DPA_gridIndex = 1 : gridPointCount          % Color each of the Grid Points according to the achieved Reliability
    
    
    
    
    gridPoint_Reliability = DPACoverageCombinationPriorityList.GPcoverageInfo{SensorCombinationIndex}(DPA_gridIndex,2);
    
    
    
    
    if ge(gridPoint_Reliability, CoverageReliabilityThreshold95) && le(gridPoint_Reliability, 1.00)
        gridPointTally95_temp(DPA_gridIndex) = 1;
    else
        gridPointTally95_temp(DPA_gridIndex) = 0;
    end
    
    
    
    if ge(gridPoint_Reliability, CoverageReliabilityThreshold50) && le(gridPoint_Reliability, 1.00)
        gridPointTally50_temp(DPA_gridIndex) = 1;
    else
        gridPointTally50_temp(DPA_gridIndex) = 0;
    end
    
    
    
    
    
    
    
    if gridPoint_Reliability ~= 0
        DPAgridPoint_ploth(DPA_gridIndex).Marker = 'o';
        DPAgridPoint_ploth(DPA_gridIndex).LineWidth = 0.5;
    end
    
    
    
    
    
    
    
    if ge(gridPoint_Reliability, 0) && lt(gridPoint_Reliability, 0.95)
        DPAgridPoint_ploth(DPA_gridIndex).Color = Rge5lt95_GPCcolor;
        DPAgridPoint_ploth(DPA_gridIndex).MarkerFaceColor = Rge5lt95_GPCcolor;
        
    elseif ge(gridPoint_Reliability, 0.95) && le(gridPoint_Reliability, 1.00)
        DPAgridPoint_ploth(DPA_gridIndex).Color = Rge95_GPcolor;
        DPAgridPoint_ploth(DPA_gridIndex).MarkerFaceColor = Rge95_GPcolor;
        
    elseif gridPoint_Reliability == 0
        DPAgridPoint_ploth(DPA_gridIndex).Marker = '.';
        DPAgridPoint_ploth(DPA_gridIndex).Color = Rlt5_GPcolor;
        DPAgridPoint_ploth(DPA_gridIndex).MarkerFaceColor = Rlt5_GPcolor;
        DPAgridPoint_ploth(DPA_gridIndex).LineWidth = 1.5;
        
        
    end     % If ends making sure grid points are colored according to the combined RC achieved
    
    
    
    
end     % For loop ends after updating the grid point information
































%%% Print the ESC and Mode selectors information
% Generating initial TextBox string
mNdiff = ECSInfo_NperLine - mod(escCountInvolved, ECSInfo_NperLine);
selectedESCsensorCombination = DPACoverageCombinationPriorityList.CombinationSensorSiteID{SensorCombinationIndex};
escCountSelectedCombination = size(selectedESCsensorCombination,1);
TxtBox_str_init_ESC = strings(escCountSelectedCombination, 1);
for escIndex = 1 : escCountSelectedCombination-1         % for each of the ESCs
    TxtBox_str_init_ESC(escIndex) = [selectedESCsensorCombination{escIndex}, ', '];
end

if escCountSelectedCombination==1
    TxtBox_str_init_ESC(escCountSelectedCombination,1) = selectedESCsensorCombination{escCountSelectedCombination};     % Collecting all the active ESC information
else
    TxtBox_str_init_ESC(escCountSelectedCombination,1) = ['and ', selectedESCsensorCombination{escCountSelectedCombination}];     % Collecting all the active ESC information
end



% TxtBox_strN_init = reshape(TxtBox_str_init_ESC, [escCountSelectedCombination, 1]);      
TxtBox_strN_init = TxtBox_str_init_ESC';

TxtBox_str_init = strings(strN_nrows+1, 1);
for strNIndex = 1 : strN_nrows+1    
    
    if strNIndex == 1        
        if sum(escTallyInvolved) == 0            
            TxtBox_str_init(strNIndex) = 'Activate ESC sensors for DPA coverage ';
        elseif sum(escTallyInvolved) == 1            
            TxtBox_str_init(strNIndex) = 'Coverage provided by the ESC sensor ';  
        else            
            TxtBox_str_init(strNIndex) = 'Coverage provided by the ESC sensor combination ';            
        end        
    else        
        TxtBox_str_init(strNIndex) = strjoin(TxtBox_strN_init(strNIndex-1, :));        
    end    
    
end     % For loop ends after printing N active ESC information per line







% Print the ESC Information: Intro and Name of the Active ESCs
hAnnotation1 = gobjects(1, strN_nrows+1);
for strNIndex = 1 : strN_nrows+1        % Plot the initial ESC Info Text    
    
    if strNIndex == 1
        hAnn1color = Statement_TB1color;
    else
        hAnn1color = Info_TB1color;
    end
    hAnnotation1(strNIndex) = annotation('textbox',TB1_Dimension_All(strNIndex, :),'String', TxtBox_str_init(strNIndex), ...
        'HorizontalAlignment', TB1_HAlignment, 'FitBoxToText','on', 'FontWeight', 'bold', 'color', hAnn1color, 'EdgeColor', 'none');    
    
end     % For loop ends after printing three (3) active ESC information per line




% Print the ESC Information: Number of Active ESCs and the Realiability of the selected Grid Point
TxBox1_str_init_end = ['Number of ESC sensors involved: ', num2str(escCountSelectedCombination)];
hAnnotation1(end+1) = annotation('textbox',TB1_Dimension_All(end, :),'String',TxBox1_str_init_end,'FitBoxToText','on', ...
    'HorizontalAlignment', TB1_HAlignment, 'color', Info_TB1color, 'EdgeColor', 'none', 'FontWeight', 'bold');








% Print the ESC Information: Number of Active ESCs and the Realiability of the selected Grid Point
hAnnotation2 = gobjects(1, 3);
TxBox2_str_init_Intro = 'ESC Sensor Combination DPA Coverage Performance:';
hAnnotation2(1) = annotation('textbox',TB2_Dimension_All(1, :),'String',TxBox2_str_init_Intro,'FitBoxToText','on', ...
    'HorizontalAlignment', TB2_HAlignment, 'color', Statement_TB2color, 'EdgeColor', 'none', 'FontWeight', 'bold');



% Print the DPA Coverage Results: Percentage of Part-1 DPA Grid Points covered with Reliability 0.95
Part1CoveragePercentage95 = sum(gridPointTally95_temp(DPA_gridInfo_Temp(:,4)==1)==1)/sum(DPA_gridInfo_Temp(:,4)==1)*100;
TxtBos_str_init_Part1CovPercentage = ['Grid Point Coverage (inside 75Km, Rel >= 0.95): ', num2str(Part1CoveragePercentage95), '% of the DPA'];
if Part1CoveragePercentage95==100
    hAnn2color = pass_TB2color;
else
    hAnn2color = fail_TB2color;
end
hAnnotation2(2) = annotation('textbox',TB2_Dimension_All(2, :),'String',TxtBos_str_init_Part1CovPercentage,'FitBoxToText','on', ...
    'HorizontalAlignment', TB2_HAlignment, 'color', hAnn2color, 'EdgeColor', 'none', 'FontWeight', 'bold');



% Print the DPA Coverage Results: Percentage of Part-2 DPA Grid Points covered with Reliability 0
Part2CoveragePercentage50 = sum(gridPointTally50_temp(DPA_gridInfo_Temp(:,4)==0)==1)/sum(DPA_gridInfo_Temp(:,4)==0)*100;
TxtBos_str_init_Part2CovPercentage = ['Grid Point Coverage (beyond 75Km, Rel >= 0): ', num2str(Part2CoveragePercentage50), '% of the DPA'];
if Part2CoveragePercentage50==100
    hAnn2color = pass_TB2color;
else
    hAnn2color = fail_TB2color;
end
hAnnotation2(3) = annotation('textbox',TB2_Dimension_All(3, :),'String',TxtBos_str_init_Part2CovPercentage,'FitBoxToText','on', ...
    'HorizontalAlignment', TB2_HAlignment, 'color', hAnn2color, 'EdgeColor', 'none', 'FontWeight', 'bold');
















% Print the sensor Indices of the first five(5) Highest Priority Sensor Combination (HPSC)
% hAnnotation3 = gobjects(1, 2);
% TxBox3_str_init_Intro = ['ESC Sensor Combination Priority Index: ', num2str(SensorCombinationIndex)];
% hAnnotation3(1) = annotation('textbox',TB3_Dimension_All(5, :),'String',TxBox3_str_init_Intro,'FitBoxToText','on', ...
%     'HorizontalAlignment', TB3_HAlignment, 'color', Statement_TB3color, 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 11);
% 
% 
% 
% selectedESCsensorCombination = DPACoverageCombinationPriorityList.CombinationSensorSiteID{SensorCombinationIndex};
% selectedCombinationsize = size(selectedESCsensorCombination,1);
% 
% if selectedCombinationsize==1
%     TxtBox3_str_init_iterationCombination = ['Sensor in the Combination: ', selectedESCsensorCombination{1}];
% else
%     TxtBox3_str_init_iterationCombination = 'Sensors in the Combination: ';
%     
%     
%     for sensorIndex = 1 : selectedCombinationsize
%         
%         if sensorIndex == size(selectedESCsensorCombination,1)
%             TxtBox3_str_init_iterationCombination = [TxtBox3_str_init_iterationCombination, 'and ' selectedESCsensorCombination{sensorIndex}];
%         else
%             TxtBox3_str_init_iterationCombination = [TxtBox3_str_init_iterationCombination, selectedESCsensorCombination{sensorIndex}, ', '];
%         end
%         
%     end
%     
% end
% 
% 
% 
% hAnn3color = HPSC_TB3color;
% hAnn3LineWidth = 1.5;
% 
% 
% hAnnotation3(2) = annotation('textbox',TB3_Dimension_All(6, :),'String',TxtBox3_str_init_iterationCombination);
% hAnnotation3(2).FitBoxToText = 'On';
% hAnnotation3(2).HorizontalAlignment = TB3_HAlignment;
% hAnnotation3(2).Color = hAnn3color;
% hAnnotation3(2).EdgeColor = hAnn3color;
% hAnnotation3(2).FontWeight = 'Bold';
% hAnnotation3(2).LineWidth = hAnn3LineWidth;










% Print the Custom Legend Information
hAnnotation4(1) = annotation('textbox', mode_dim1_all(1,:), 'String', 'O  95%-100% Reliability', 'FitBoxToText', 'on', ...
    'HorizontalAlignment', MS_HAlignment, 'color', [0, 0.5, 0], 'FontWeight', 'bold', 'LineStyle', 'none');
hAnnotation4(2) = annotation('textbox', mode_dim1_all(2,:), 'String', 'O  50%-95% Reliability', 'FitBoxToText', 'on', ...
    'HorizontalAlignment', MS_HAlignment, 'color', 'b', 'FontWeight', 'bold', 'LineStyle', 'none');
hAnnotation4(3) = annotation('textbox', mode_dim1_all(3,:), 'String', 'O  <50% Reliability', 'FitBoxToText', 'on', ...
    'HorizontalAlignment', MS_HAlignment, 'color', 'r', 'FontWeight', 'bold', 'LineStyle', 'none');
customLegend_BorderDim = [mode_dim1_all(end,1:3), mode_dim1_all(1,2)-mode_dim1_all(end,2)+mode_dim1_all(end,4)];
hAnnotation4(4) = annotation('textbox', customLegend_BorderDim, 'String', ' ', 'FitBoxToText', 'off', ...
    'HorizontalAlignment', MS_HAlignment, 'EdgeColor', 'k', 'LineWidth', 1);










% Print the title and axis labels of the ESC Logic System visualization
if ischar(DPA_IDstr)    
    caption = sprintf('ESC Logic System: Analysis of ESC sensor combination selection priority:    DPA Zone: %s,  DPA ID: %s', DPA_Zone, DPA_IDstr);
else
    caption = sprintf('ESC Logic System: Analysis of ESC sensor combination selection priority:    DPA Zone: %s,  DPA ID: %d', DPA_Zone, DPA_IDstr);
end
title(caption);       
xlabel('Longitude', 'FontWeight', 'bold'); ylabel('Latitude', 'FontWeight', 'bold');
ax = gca; ax.FontWeight = 'bold';
























xx=1;
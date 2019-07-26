%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_GenerateESCsensorInformation.m
%%%%
%%%%        - Function that Generates and Saves ESC sensor Info associated with the given DPA
%%%%            - ESC sensor Location Information
%%%%            - ESC sensor Antenna Information
%%%%
%%%%        - Generates and saves ESC sensor Information Table
%%%%            - DPAName_ESCsensorInfo.xlsx
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function ESCsensorInfo = Function_GenerateESCsensorInformation(DPA_Zone, DPA_ID)



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






DPA_Name = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
DPArelated_dataSaveFolder = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Original DPA and ESC Data';
DPAspecific_dataSaveFolder = strcat(DPArelated_dataSaveFolder, '\', DPA_Name);
if exist(DPAspecific_dataSaveFolder, 'dir') ~= 7
    mkdir(DPAspecific_dataSaveFolder)
end













xx=1;   % Completes Input initialization














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Gather and Save the ESC sensor Information for the given DPA
%%%%    ESC_Info: for each of the N ESC sensors associated with the DPA, Struct with fields:
%%%%            - Longitude: [Longitude] (Nx1 Double)
%%%%            - Latitude: [Latitude] (Nx1 Double)
%%%%            - Antenna Height: [Height] (Nx1 Double)
%%%%            - Antenna Azimuth: [AzimuthAngle] (Nx1 Double)
%%%%            - Site ID: [SiteID] (Nx1 String)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% Load ESC sensor Information for the given DPA
% DPAESCLocation_folderName = 'C:\Users\MunawwarSohul\Dropbox (Federated Wireless)\007 Simulations\003 Related FW Documents\Data for Analysis';
% ESC_dataFolderName = strcat(DPAESCLocation_folderName, '\007 FW DPA ESC Information');
% % ESC_Data_fileName = 'FW_DPA_ESC_Locations_POR_02_19_19.xlsx';
% ESC_Data_fileName = 'FW_DPA_ESC_Locations_POR_02_19_19_updated.xlsx';
% ESC_Data_fullFileName = fullfile(ESC_dataFolderName, ESC_Data_fileName);
% 
% SheetName = 'Sheet1';
% ESC_Data = readtable(ESC_Data_fullFileName, 'sheet', SheetName);
% escInfoIndex = find(strcmp(ESC_Data.DPAID, DPA_Name) == 1);
% 
% ESCsensorInfo.SiteID = ESC_Data.SiteID(escInfoIndex);
% ESCsensorInfo.Longitude = ESC_Data.Longitude(escInfoIndex);
% ESCsensorInfo.Latitude = ESC_Data.Latitude(escInfoIndex);
% ESCsensorInfo.AntennaHeight = ESC_Data.ESCAGLInFeet(escInfoIndex) * 0.3048;
% % ESCInfo.AntennaAzimuth = ESC_Data.Azimuth(escInfoIndex);
% ESCsensorInfo.AntennaAzimuth = ESC_Data.Azimuth_updated(escInfoIndex); %#ok<STRNU>
% 
% escCount = length(ESCsensorInfo.SiteID);
% ESCInfoDPAname = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
% escInfoIndex = find(strcmp(ESC_Data.DPAID, ESCInfoDPAname) == 1);
% 
% ESCsensorInfo.SensorCategory = zeros(size(ESCsensorInfo.SiteID));
% ESCsensorInfo.SensorCategory(escInfoIndex) = ones(size(escInfoIndex));
% 
% 
% % Save ESCsensorInfo
% ESCsensorInfo_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_ESCsensorInfo.mat');
% ESCsensorInfo_fullFileName = fullfile(DPAspecific_dataSaveFolder, ESCsensorInfo_fileName);
% if exist(ESCsensorInfo_fullFileName, 'file') == 2
%     delete(ESCsensorInfo_fullFileName)
% end
% save(ESCsensorInfo_fullFileName, 'ESCsensorInfo')
% clear ESCsensorInfo_fileName ESCsensorInfo_fullFileName










%%% load ESC information associated with the ESC logic System analysis
ELS_ESCInfo_folderName = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Additional Data for Analysis';
if strcmp(DPA_Zone, 'east')
    ELS_ESCInfo_fullFileName = strcat(ELS_ESCInfo_folderName, '\FW_ELS_ESC_Information_06_03_19_EAST.xlsx');
elseif strcmp(DPA_Zone, 'west')
    ELS_ESCInfo_fullFileName = strcat(ELS_ESCInfo_folderName, '\FW_ELS_ESC_Information_06_03_19_WEST.xlsx');
end

ESCInfoDPAname = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
ESCsensor_Data = readtable(ELS_ESCInfo_fullFileName, 'sheet', ESCInfoDPAname);


% Gather the ESC Sensor Information
ESCsensorInfo.sensorIndex = (1:size(ESCsensor_Data,1))';
ESCsensorInfo.SiteID = ESCsensor_Data.SiteID;
ESCsensorInfo.Longitude = ESCsensor_Data.Longitude;
ESCsensorInfo.Latitude = ESCsensor_Data.Latitude;
ESCsensorInfo.AntennaHeight = ESCsensor_Data.HeightAGL * 0.3048;
ESCsensorInfo.Azimuth = ESCsensor_Data.Azimuth;
ESCsensorInfo.Elevation = zeros(size(ESCsensor_Data,1),1);

ESCInfoDPAname = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
escInfoIndex = find(strcmp(ESCsensor_Data.DPAID, ESCInfoDPAname) == 1);

ESCsensorInfo.SensorCategory = zeros(size(ESCsensorInfo.SiteID));
ESCsensorInfo.SensorCategory(escInfoIndex) = ones(size(escInfoIndex));



% Save ESCsensorInfo
ESCsensorInfo_fileName = strcat(DPA_Zone, '_', num2str(DPA_IDstr), '_ESCsensorInfo.mat');
ESCsensorInfo_fullFileName = fullfile(DPAspecific_dataSaveFolder, ESCsensorInfo_fileName);
if exist(ESCsensorInfo_fullFileName, 'file') == 2
    delete(ESCsensorInfo_fullFileName)
end
save(ESCsensorInfo_fullFileName, 'ESCsensorInfo')
clear ESCsensorInfo_fileName ESCsensorInfo_fullFileName



% Create the ESC Sensor Info table to be saved as .xlsx file
ESCsensorInfoTable = table(ESCsensorInfo.sensorIndex, ESCsensorInfo.Longitude, ESCsensorInfo.Latitude, ...
    ESCsensorInfo.AntennaHeight, ESCsensorInfo.Azimuth, ESCsensorInfo.Elevation);

ESCsensorInfoTable.Properties.VariableNames{1} = 'ESCsensorIndex';
ESCsensorInfoTable.Properties.VariableNames{2} = 'ESCsensorLon';
ESCsensorInfoTable.Properties.VariableNames{3} = 'ESCsensorLat';
ESCsensorInfoTable.Properties.VariableNames{4} = 'ESCsensorHeight';
ESCsensorInfoTable.Properties.VariableNames{5} = 'ESCsensorAzimuth';
ESCsensorInfoTable.Properties.VariableNames{6} = 'ESCsensorElevation';


% Save the ESC sensor Information (.xlsx file)
DPAfolderName = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
ELS_dataSaveFolder = strcat('C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Data_20190531\', DPAfolderName);
ESCsensorInfo_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_ESCsensorInfo.xlsx');
ESCsensorInfo_fullFileName = fullfile(ELS_dataSaveFolder, ESCsensorInfo_filename);
if exist(ESCsensorInfo_fullFileName, 'file') == 2
    delete(ESCsensorInfo_fullFileName)
end

writetable(ESCsensorInfoTable, ESCsensorInfo_fullFileName, 'Sheet', 1, 'Range','A1')
clear ESCsensorInfo_filename ESCsensorInfo_fullFileName




















xx=1;       fprintf('\nESC Logic System Analysis: Generation of ESC Sensor Information :: Complete\n')











end
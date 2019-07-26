%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Function_AntennaPatternProcessing.m
%%%%
%%%%        - Function that processes and saves the Antenna Gain Pattern
%%%%            - Resolution of the gain pattern is changed from 3degree to 1 degree
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function Function_GenerateAntennaPattern_MultipleESC(DPA_Zone, DPA_ID, ESCsensorInfo)



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







DPAfolderName = strcat(DPA_Zone, '_dpa_', num2str(DPA_IDstr));
ELS_dataSaveFolder = strcat('C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Data_20190531\', DPAfolderName);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESC Sensor Antenna Gain Pattern Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELS_AntennaPatternInfo_folderName = 'C:\Users\czhu\Matlab simulation zch\Azimuth Optimization\Additional Data for Analysis';
ELS_AntennaPatternInfo_fullFileName = strcat(ELS_AntennaPatternInfo_folderName, '\FW_ESC_AntennaGainPatternInfo.xlsx');

AntennaPattern_Data = readtable(ELS_AntennaPatternInfo_fullFileName, 'sheet', 'Sheet1');
AntennaGainValue = table2array(AntennaPattern_Data);


escCount = size(ESCsensorInfo.sensorIndex,1);
escAntennaPattern = AntennaGainValue(:,1);
for escIndex = 1 : escCount     % For each ESC sensors    
    
    % Get the ESC Antenna Pattern
    AntennaRotation = mod(roundn((AntennaGainValue(:,1) + ESCsensorInfo.Azimuth(escIndex)), 1), 360);
    AntennaRotationIndex = find(AntennaGainValue(:,1)==AntennaRotation(1)) - 1;
    escAntennaPattern(:,escIndex+1) = circshift(AntennaGainValue(:,2), AntennaRotationIndex);    
    
end





%%% Create the Antenna Pattern Info table to be saved as .xlsx file
AntennaPatternInfoTable = array2table(escAntennaPattern);
AntennaPatternInfoTable.Properties.VariableNames{1} = 'AzimuthAngle';

for escIndex = 1 : escCount
    AntennaPatternInfoTable.Properties.VariableNames{escIndex+1} = strcat('ESCsensor', num2str(escIndex), '_Gain');
end

% Save the Antenna Pattern Information (.xlsx file)
AntennaPatternInfo_filename = strcat(DPA_Zone, '_', num2str(DPA_ID), '_AntennaPatternInfo.xlsx');
AntennaPatternInfo_fullFileName = fullfile(ELS_dataSaveFolder, AntennaPatternInfo_filename);

writetable(AntennaPatternInfoTable, AntennaPatternInfo_fullFileName, 'Sheet', 1, 'Range','A1')








end
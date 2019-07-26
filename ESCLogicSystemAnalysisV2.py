#
# ESC Logic System Analysis
#

import numpy as np
import pandas as pd
# from pandas import ExcelWriter
# from pandas import ExcelFile
# import copy
# import sys
# import math
from src.harness.reference_models.propagation.wf_itm import CalcItmPropagationLoss
from src.harness.reference_models.geo.vincenty import GeodesicDistanceBearing
# from src.harness.reference_models.geo.vincenty import GeodesicPoint


# Simple round up function
# def roundup(x):
#     return int(math.ceil(x / 10.0)) * 10


# Setup the file names for each of the DPAs
DPA_zone = 'east'
DPA_ID = 1

DPA_folderName = DPA_zone+'_dpa_'+str(DPA_ID)
DPA_Name = DPA_zone+'_'+str(DPA_ID)
DPAGridPintInfo_fileName = 'Data_20190531/'+DPA_folderName+'/'+DPA_Name+'_DPAGridPointInfo.xlsx'
ESCsensorInfo_fileName = 'Data_20190531/'+DPA_folderName+'/'+DPA_Name+'_ESCsensorInfo.xlsx'
AntennaPatternInfo_fileName = 'Data_20190531/'+DPA_folderName+'/'+DPA_Name+'_AntennaPatternInfo.xlsx'
GridPointDetectionInfo_fileName = 'Data_20190531/'+DPA_folderName+'/'+DPA_Name+'_DPAGridPointDetectionInfo.xlsx'
GridPointCoverageInfo_fileName = 'Data_20190531/'+DPA_folderName+'/'+DPA_Name+'_DPAGridPointCoverageInfo.xlsx'


# Get the DPA Grid Point Information
DPAGridPintInfo_DatafromFile = pd.ExcelFile(DPAGridPintInfo_fileName)
DPAGridPintInfo_DataFrame = DPAGridPintInfo_DatafromFile.parse('Sheet1')

GridPointinfo_VariableNames = DPAGridPintInfo_DataFrame.keys()
GridPointIndex = DPAGridPintInfo_DataFrame[GridPointinfo_VariableNames[0]]
GridPointLon = DPAGridPintInfo_DataFrame[GridPointinfo_VariableNames[1]]
GridPointLat = DPAGridPintInfo_DataFrame[GridPointinfo_VariableNames[2]]
GridPointInnerBoundaryIndex = DPAGridPintInfo_DataFrame[GridPointinfo_VariableNames[3]]

GridPointInnerBoundaryIndexArray = np.asarray(GridPointInnerBoundaryIndex)

# Get the ESC sensor Information
ESCsensorInfo_DatafromFile = pd.ExcelFile(ESCsensorInfo_fileName)
ESCsensorInfo_DataFrame = ESCsensorInfo_DatafromFile.parse('Sheet1')

ESCsensorInfo_VariableNames = ESCsensorInfo_DataFrame.keys()
ESCsensorIndex = ESCsensorInfo_DataFrame[ESCsensorInfo_VariableNames[0]]
ESCsensorLon = ESCsensorInfo_DataFrame[ESCsensorInfo_VariableNames[1]]
ESCsensorLat = ESCsensorInfo_DataFrame[ESCsensorInfo_VariableNames[2]]
ESCsensorHeight = ESCsensorInfo_DataFrame[ESCsensorInfo_VariableNames[3]]
ESCsensorAzimuth = ESCsensorInfo_DataFrame[ESCsensorInfo_VariableNames[4]]
ESCsensorElevation = ESCsensorInfo_DataFrame[ESCsensorInfo_VariableNames[5]]


# Get the Antenna Pattern Information
AntennaPatternInfo_DatafromFile = pd.ExcelFile(AntennaPatternInfo_fileName)
AntennaPatternInfo_DataFrame = AntennaPatternInfo_DatafromFile.parse('Sheet1')

AntennaPatternInfo_VariableNames = AntennaPatternInfo_DataFrame.keys()
AntennaPattern_AzimuthAngle = AntennaPatternInfo_DataFrame[AntennaPatternInfo_VariableNames[0]]


# For each of the ESC sensors generate the GridPoint Detection and Coverage information
radarTxEIRP = 121
escTargetRSSI = -89
ReliabilityList = np.arange(0.5, 1.0, 0.01)

InsideInnerBoundary_ReliabilityThreshold = 0.95
OutsideInnerBoundary_ReliabilityThreshold = 0.50

GridPointDetectionInfo = np.empty([GridPointIndex.shape[0], ESCsensorIndex.shape[0]], dtype=float)
for gpIndex in GridPointIndex:

    for escIndex in ESCsensorIndex:

        # Find the ESC Rx antenna gain towards the DPA grid point
        DistanceAzimuthInfo = GeodesicDistanceBearing(ESCsensorLat[escIndex-1], ESCsensorLon[escIndex-1],
                                                      GridPointLat[gpIndex-1], GridPointLon[gpIndex-1],
                                                      accuracy=1.0E-12)
        azimuthAngle = DistanceAzimuthInfo[1]
        azimuthAngleRounded = round(azimuthAngle)
        if azimuthAngleRounded == 360:
            azimuthAngleRounded = 0

        escAntennaGainAzimuthIndex, = np.where(AntennaPattern_AzimuthAngle == azimuthAngleRounded)
        escAntennaGainAzimuthListArray = np.asarray(AntennaPatternInfo_DataFrame)
        escAntennaGainAzimuth = escAntennaGainAzimuthListArray[escAntennaGainAzimuthIndex, escIndex]

        # Get the PathLoss values for the current ESC sensor and grid point for all the reliability values
        PLcbsd = CalcItmPropagationLoss(GridPointLat[gpIndex-1], GridPointLon[gpIndex-1], 50,
                                        ESCsensorLat[escIndex-1], ESCsensorLon[escIndex-1], ESCsensorHeight[escIndex-1],
                                        False, ReliabilityList,
                                        freq_mhz=3625.,
                                        its_elev=None,
                                        is_height_cbsd_amsl=False,
                                        return_internals=False)
        PathLoss_RefCode = PLcbsd.db_loss

        # Calculate the received signal power at the ESC sensor
        escRcvdRSSI = [radarTxEIRP - PathLossValue + escAntennaGainAzimuth[0] for PathLossValue in PathLoss_RefCode]

        # Find the DPA grid point detection and coverage information
        escRcvdRSSIarray = np.asarray(escRcvdRSSI)
        if any(escRcvdRSSIarray > escTargetRSSI):

            # Find the maximum reliability with which the current ESC sensor detects the current gridpoint
            maxReliabilityIndex = max(np.argwhere(escRcvdRSSIarray > escTargetRSSI))

            # Record the Grid Point Detection Information
            GridPointDetectionInfo[gpIndex-1, escIndex-1] = ReliabilityList[maxReliabilityIndex]

        xx = 1

    xx1 = 1

# Write the DPA Grid Point Detection and Coverage Information to Excel file
GridPointDetectionInfo_DataFrame = pd.DataFrame(GridPointDetectionInfo)
GridPointDetectionInfo_ExcelHandle = pd.ExcelWriter(GridPointDetectionInfo_fileName, engine='xlsxwriter')
GridPointDetectionInfo_DataFrame.to_excel(GridPointDetectionInfo_ExcelHandle, sheet_name='Sheet1')
GridPointDetectionInfo_ExcelHandle.save()

xx = 1

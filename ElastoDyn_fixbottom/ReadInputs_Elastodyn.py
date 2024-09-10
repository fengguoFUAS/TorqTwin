#　This function ......................................
# Author: Dr.Ing. Feng Guo (郭峰)
# Contacts: fengguokm@outlook.com
# licensed under the Apache License 2.0.

import numpy as np

class InitInputData:
    def __init__(self,DOF_Flag,IniCondition,TurbConfig,Blade,Tower):
        self.DOF_Flag = DOF_Flag  # for saving DOF flags
        self.IniCondition = IniCondition  # for initial conditions
        self.TurbConfig = TurbConfig  # for turbine configurations, I also put the mass and inertia, drivetrain info into TurbConfig
        self.Blade = Blade # for blade input data
        self.Tower = Tower  # for tower input data
        self.TurbConfig['TwrNodes'] = Tower['TwrNodes']
        self.TurbConfig['BldNodes'] = Blade['BldNodes']
        self.TurbConfig['RotSpeed0'] = IniCondition['RotSpeed']


def ReadInputs_Elastodyn(ElastoDynInputFileName):
    DOF_Flag      = {}
    IniCondition  = {}
    TurbConfig    = {}
    Blade         = {}
    Tower         = {}

    DOF_string    = {'FlapDOF1','FlapDOF2','EdgeDOF','TeetDOF','DrTrDOF','GenDOF','YawDOF','TwFADOF1','TwFADOF2',
                    'TwSSDOF1','TwSSDOF2','PtfmSgDOF','PtfmSwDOF','PtfmHvDOF','PtfmRDOF','PtfmPDOF','PtfmYDOF'}
    IniCondi_string = {'OoPDefl','IPDefl','BlPitch(1)','BlPitch(2)','BlPitch(3)','TeetDefl','Azimuth','RotSpeed',
                       'NacYaw','TTDspFA','TTDspSS','PtfmSurge','PtfmSway','PtfmHeave','PtfmRoll','PtfmPitch','PtfmYaw'}
    TurbConfig_string ={'NumBl','TipRad','HubRad','PreCone(1)','PreCone(2)','PreCone(3)','HubCM','UndSling','Delta3',
                        'AzimB1Up',	'OverHang','ShftGagL','ShftTilt','NacCMxn','NacCMyn','NacCMzn','NcIMUxn','NcIMUyn',
                        'NcIMUzn','Twr2Shft','TowerHt','TowerBsHt',	'PtfmCMxt',	'PtfmCMyt',	'PtfmCMzt',	'PtfmRefzt',
                        'TipMass(1)','TipMass(2)','TipMass(3)','HubMass','HubIner','GenIner','NacMass',	'NacYIner',
                        'YawBrMass','PtfmMass','PtfmRIner','PtfmPIner','PtfmYIner','GBoxEff','GBRatio','DTTorSpr','DTTorDmp'}

    with open(ElastoDynInputFileName, 'r') as file:

        for current_line in file:

            if 'SumPrint' in current_line:  # the lines in the output section are ignored
                break

            for DOF_tosearch in DOF_string:
                if DOF_tosearch in current_line:
                    DOF_Flag[DOF_tosearch] = current_line.split(DOF_tosearch)[0].strip()

            for IniCondi_tosearch in IniCondi_string:
                if IniCondi_tosearch in current_line and 'AzimB1Up' not in current_line:# the description includes the variable name
                    IniCondition[IniCondi_tosearch] = float(current_line.split(IniCondi_tosearch)[0].strip())

            for TurbConfig_tosearch in TurbConfig_string:
                if TurbConfig_tosearch in current_line:
                    TurbConfig[TurbConfig_tosearch] = float(current_line.split(TurbConfig_tosearch)[0].strip())

            if 'BldNodes' in current_line in current_line:
                Blade['BldNodes'] = int(current_line.split('BldNodes')[0].strip())
            if 'BldFile1' in current_line:  # some input files uses BldFile1 but some uses BldFile(1) as the identifier
                Blade['BldFile'] = current_line.split('BldFile1')[0].strip().replace('"', "")   # the identical blades are assumed
            if 'BldFile(1)' in current_line:
                Blade['BldFile'] = current_line.split('BldFile(1)')[0].strip().replace('"', "")   # the identical blades are assumed
            if 'TwrNodes' in current_line in current_line:  # the description includes the variable name
                Tower['TwrNodes'] = int(current_line.split('TwrNodes')[0].strip())
            if 'TwrFile' in current_line:
                Tower['TwrFile'] = current_line.split('TwrFile')[0].strip().replace('"', "")

    InputData=InitInputData(DOF_Flag,IniCondition,TurbConfig,Blade,Tower)
    return InputData


def ReadInputs_Blade(BladeInputFileName):  #BLADE ADJUSTMENT FACTORS are ingored
    BladeData = {}
    BladeData['BldFl1Sh'] = np.zeros(6, dtype=float)
    BladeData['BldFl2Sh'] = np.zeros(6, dtype=float)
    BladeData['BldEdgSh'] = np.zeros(6, dtype=float)

    with open(BladeInputFileName, 'r') as file:

        for current_line in file:
            if 'NBlInpSt' in current_line:
                BladeData['NBlInpSt'] = int(current_line.split('NBlInpSt')[0].strip())
                BladeData['BlFract']   = np.zeros(BladeData['NBlInpSt'], dtype=float)
                BladeData['PitchAxis'] = np.zeros(BladeData['NBlInpSt'], dtype=float)
                BladeData['StrcTwst']  = np.zeros(BladeData['NBlInpSt'], dtype=float)
                BladeData['BMassDen']  = np.zeros(BladeData['NBlInpSt'], dtype=float)
                BladeData['FlpStff']   = np.zeros(BladeData['NBlInpSt'], dtype=float)
                BladeData['EdgStff']   = np.zeros(BladeData['NBlInpSt'], dtype=float)

            if 'BldFlDmp1' in current_line:
                BladeData['BldFlDmp1'] = float(current_line.split('BldFlDmp1')[0].strip())
            if 'BldFlDmp2' in current_line:
                BladeData['BldFlDmp2'] = float(current_line.split('BldFlDmp2')[0].strip())
            if 'BldEdDmp1' in current_line:
                BladeData['BldEdDmp1'] = float(current_line.split('BldEdDmp1')[0].strip())

            if 'DISTRIBUTED BLADE PROPERTIES' in current_line:
                current_line = file.readline().strip()
                current_line = file.readline().strip()

                for counter in range(BladeData['NBlInpSt']):
                    current_line = file.readline().strip().split()

                    Number_array = [float(num) for num in current_line]
                    BladeData['BlFract'][counter]          = Number_array[0]
                    BladeData['PitchAxis'][counter]        = Number_array[1]
                    BladeData['StrcTwst'][counter]         = Number_array[2]
                    BladeData['BMassDen'][counter]         = Number_array[3]
                    BladeData['FlpStff'][counter]          = Number_array[4]
                    BladeData['EdgStff'][counter]          = Number_array[5]

            if 'BLADE MODE SHAPES' in current_line:
                for counter in range(5):
                    current_line = file.readline().strip()
                    BladeData['BldFl1Sh'][counter+1] = float(current_line.split('BldFl1Sh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)
                for counter in range(5):
                    current_line = file.readline().strip()
                    BladeData['BldFl2Sh'][counter+1] = float(current_line.split('BldFl2Sh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)
                for counter in range(5):
                    current_line = file.readline().strip()
                    BladeData['BldEdgSh'][counter+1] = float(current_line.split('BldEdgSh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)

    return BladeData

def ReadInputs_Tower(TowerInputFileName):    #TOWER ADJUSTMUNT FACTORS are ingored
    TowerData = {}
    TowerData['TwFAM1Sh'] = np.zeros(6, dtype=float)
    TowerData['TwFAM2Sh'] = np.zeros(6, dtype=float)
    TowerData['TwSSM1Sh'] = np.zeros(6, dtype=float)
    TowerData['TwSSM2Sh'] = np.zeros(6, dtype=float)

    with open(TowerInputFileName, 'r') as file:

        for current_line in file:
            if 'NTwInpSt' in current_line:
                TowerData['NTwInpSt']  = int(current_line.split('NTwInpSt')[0].strip())
                TowerData['HtFract']   = np.zeros(TowerData['NTwInpSt'], dtype=float)
                TowerData['TMassDen']  = np.zeros(TowerData['NTwInpSt'], dtype=float)
                TowerData['TwFAStif']  = np.zeros(TowerData['NTwInpSt'], dtype=float)
                TowerData['TwSSStif']  = np.zeros(TowerData['NTwInpSt'], dtype=float)


            if 'TwrFADmp(1)' in current_line:
                TowerData['TwrFADmp(1)'] = float(current_line.split('TwrFADmp(1)')[0].strip())
            if 'TwrFADmp(2)' in current_line:
                TowerData['TwrFADmp(2)'] = float(current_line.split('TwrFADmp(2)')[0].strip())
            if 'TwrSSDmp(1)' in current_line:
                TowerData['TwrSSDmp(1)'] = float(current_line.split('TwrSSDmp(1)')[0].strip())
            if 'TwrSSDmp(2)' in current_line:
                TowerData['TwrSSDmp(2)'] = float(current_line.split('TwrSSDmp(2)')[0].strip())

            if 'DISTRIBUTED TOWER PROPERTIES' in current_line:
                current_line = file.readline().strip()
                current_line = file.readline().strip()

                for counter in range(TowerData['NTwInpSt']):
                    current_line = file.readline().strip().split()

                    Number_array = [float(num) for num in current_line]
                    TowerData['HtFract'][counter]  = Number_array[0]
                    TowerData['TMassDen'][counter] = Number_array[1]
                    TowerData['TwFAStif'][counter] = Number_array[2]
                    TowerData['TwSSStif'] [counter]= Number_array[3]

            if 'TOWER FORE-AFT MODE SHAPES' in current_line:
                for counter in range(5):
                    current_line = file.readline().strip()
                    TowerData['TwFAM1Sh'][counter+1] = float(current_line.split('TwFAM1Sh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)
                for counter in range(5):
                    current_line = file.readline().strip()
                    TowerData['TwFAM2Sh'][counter+1] = float(current_line.split('TwFAM2Sh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)
            if 'TOWER SIDE-TO-SIDE MODE SHAPES' in current_line:
                for counter in range(5):
                    current_line = file.readline().strip()
                    TowerData['TwSSM1Sh'][counter+1] = float(current_line.split('TwSSM1Sh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)
                for counter in range(5):
                    current_line = file.readline().strip()
                    TowerData['TwSSM2Sh'][counter+1] = float(current_line.split('TwSSM2Sh(' + str(counter+2)+ ')')[0].strip())  #　coefficient of x^(counter+2)

    return TowerData
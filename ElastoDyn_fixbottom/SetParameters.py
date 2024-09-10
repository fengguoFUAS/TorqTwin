import numpy as np

deg2rad = np.deg2rad(1)

def Shape(Fract, FlexL, ModShpAry, Deriv):
    # Initialize Swtch
    Swtch = [0, 0, 0]

    # Check conditions for Deriv and Fract
    if not (0 <= Deriv <= 2):
        raise ValueError('Function SHP input Deriv is invalid. Deriv must be 0, 1, or 2.')
    elif not (0.0 <= Fract <= 1.0):
        raise ValueError('Function SHP input Fract does not meet the condition 0 <= Fract <= 1.')

    # Set Swtch based on Deriv value
    Swtch[Deriv] = 1

    # Initialize ShapeFunc
    ShapeFunc = 0.0

    # Compute ShapeFunc based on conditions
    for I in range(1,len(ModShpAry)):  # Assuming ModShpAry is defined elsewhere
        J = I + 1
        CoefTmp = Swtch[0] + Swtch[1] * J + Swtch[2] * I * J

        if J == 2 and Deriv == 2:
            ShapeFunc = ModShpAry[I ] * CoefTmp / (FlexL ** Deriv)
        else:
            ShapeFunc += ModShpAry[I] * CoefTmp * (Fract ** (J - Deriv)) / (FlexL ** Deriv)

    return ShapeFunc

class SetParameters:
    def __init__(self,Platform,Tower,Nacelle,Hub,Blade,DriveTrain,DOFs):
        self.Platform   = Platform
        self.Tower      = Tower
        self.Nacelle    = Nacelle
        self.Hub        = Hub
        self.Blade      = Blade
        self.DriveTrain = DriveTrain
        self.DOFs       = DOFs
        self.Gravity    = 9.81


class SetPlatformPara:
    def __init__(self,TurbConfig):
        #  distance from the ground level [onshore], MSL [offshore wind] to the platform Center of Mass [meters]
        self.PtfmCM_xyz = np.array([TurbConfig['PtfmCMxt'], TurbConfig['PtfmCMyt'], TurbConfig['PtfmCMzt']])
        self.PtfmRefzt  = TurbConfig['PtfmRefzt']
        self.PtfmMass   = TurbConfig['PtfmMass']
        self.PtfmRIner  = TurbConfig['PtfmRIner']
        self.PtfmPIner  = TurbConfig['PtfmPIner']
        self.PtfmYIner  = TurbConfig['PtfmYIner']

class SetNacellePara:
    def __init__(self, TurbConfig):

        self.NacMass  = TurbConfig['NacMass']
        self.NacYIner = TurbConfig['NacYIner']
        self.YawBrMass = TurbConfig['YawBrMass']

        # distance from the tower-top to the nacelle CM (meters)
        self.NacCM_xyz = np.array([TurbConfig['NacCMxn'], TurbConfig['NacCMyn'], TurbConfig['NacCMzn']])

class SetDriveTrainPara:
    def __init__(self, TurbConfig):
        #  Downwind distance from the ground level [onshore], MSL [offshore wind] to the platform Center of Mass [meters]
        self.Twr2Shft = TurbConfig['Twr2Shft']
        self.ShftTilt = TurbConfig['ShftTilt']
        self.GBoxEff  = TurbConfig['GBoxEff']
        self.GBRatio  = TurbConfig['GBRatio']
        self.GeneIner = TurbConfig['GenIner']
        self.DTTorSpr = TurbConfig['DTTorSpr']
        self.DTTorDmp = TurbConfig['DTTorDmp']
        self.RotSpeed0 = TurbConfig['RotSpeed0']/60*2*np.pi

class SetHubPara:
    def __init__(self, TurbConfig):
        self.HubMass = TurbConfig['HubMass']
        self.HubIner = TurbConfig['HubIner']
        self.HubCM   = TurbConfig['HubCM']
        self.HubRad  = TurbConfig['HubRad']
        self.OverHang = TurbConfig['OverHang']
        self.HubHt   = TurbConfig['TowerHt']+TurbConfig['Twr2Shft']+TurbConfig['OverHang']*np.sin(TurbConfig['ShftTilt']*deg2rad)


class SetBladePara:
    def __init__(self, TurbConfig,BladeData):
        self.NumBl    = TurbConfig['NumBl']
        self.TipRad   = TurbConfig['TipRad']
        self.BldFlexL = TurbConfig['TipRad']-TurbConfig['HubRad']
        self.PreCone  = TurbConfig['PreCone(1)']
        self.AzimB1Up = TurbConfig['AzimB1Up']
        self.BldNodes = TurbConfig['BldNodes']

        self.DRNodes     = np.zeros((self.BldNodes, 1))
        self.DRNodes[:]  = self.BldFlexL / self.BldNodes
        self.RNodes      = np.zeros((self.BldNodes, 1))
        self.RNodes[0]   = 0.5 * self.DRNodes[0]
        for J in range(1,self.BldNodes): self.RNodes[J] = self.RNodes[J - 1] + 0.5 * (self.DRNodes[J] + self.DRNodes[J - 1])
        self.RNodesNorm  = self.RNodes/self.BldFlexL
        self.TipNode     = self.BldNodes+1

        # define blade properties, identical blades are assumed
        self.Theta_s       = np.zeros((self.TipNode+1,1))
        self.Theta_s[0]    = BladeData['StrcTwst'][0]
        self.Theta_s[-1]   = BladeData['StrcTwst'][-1]
        self.Theta_s[1:-1] = np.interp(self.RNodesNorm, BladeData['BlFract'], BladeData['StrcTwst'])
        self.PitchAxis     = np.interp(self.RNodesNorm, BladeData['BlFract'], BladeData['PitchAxis'])
        self.MassB         = np.interp(self.RNodesNorm, BladeData['BlFract'], BladeData['BMassDen'])
        self.StiffBF       = np.interp(self.RNodesNorm, BladeData['BlFract'], BladeData['FlpStff'])
        self.StiffBE       = np.interp(self.RNodesNorm, BladeData['BlFract'], BladeData['EdgStff'])
        self.BldEdgSh      = BladeData['BldEdgSh']
        self.BldFl1Sh      = BladeData['BldFl1Sh']
        self.BldFl2Sh      = BladeData['BldFl2Sh']

        # calculate the mass of each blade element
        self.BElmntMass    = np.zeros((self.BldNodes, 1))
        for j in range(self.BldNodes) : self.BElmntMass[j]=self.MassB[j]*self.DRNodes[j]

        # calculate the first moment of blade section about the rotating axis above element j
        # only needed if the centrifugal stiffing effect is considered
        # FMomAbvNd = np.zeros((self.BldNodes, 1))
        # for j in range(self.BldNodes-1,-1,-1):
        #     FMomAbvNd[j] = (0.5 * self.BElmntMass[j]) * (TurbConfig['HubRad'] + self.RNodes[j] + 0.5 * self.DRNodes[j])
        #     if j == self.BldNodes-1:
        #         FMomAbvNd[j] = FMomAbvNd[j] # the tip mass is ignored
        #     else:
        #         FMomAbvNd[j] = FMomAbvNd[j]+FMomAbvNd[j+1]+(0.5 * self.BElmntMass[j+1]) * (TurbConfig['HubRad'] + self.RNodes[j+1] + 0.5 * self.DRNodes[j+1])

        # allocation
        self.GenMassBF      = np.zeros((2,1)) # generalized mass for blade flap-wise motion
        self.GenMassBE      = np.zeros((1,1)) # generalized mass for blade edge-wise motion
        self.GenStiffBF     = np.zeros((2,2)) # generalized stiffness for blade flap-wise motion
        self.GenStiffBE     = np.zeros((1, 1))  # generalized stiffness for blade edge-wise motion
        self.GenDampBF      = np.zeros((2, 2))  # generalized stiffness for blade flap-wise motion
        self.GenDampBE      = np.zeros((1, 1))  # generalized stiffness for blade edge-wise motion
        self.FreqBF         = np.zeros((2,1))   # natural frequencies of blade flap-wise mode without centrifugal effect
        self.FreqBE         = np.zeros((1, 1))  # natural frequencies of blade edge-wise mode without centrifugal effect

        #self.CenStiffBF     = np.zeros((2, 2))  # centrifugal stiffness for blade flap-wise motion
        #self.CenStiffBE     = np.zeros((1, 1))  # centrifugal stiffness for blade edge-wise motion
        self.TwistedSF      = np.zeros((2,3,self.TipNode+1,3)) # twisted shape function, (sin+cos,three DOF, BldNodes+2,0+1+2 derivatives)
        self.AxRedBldSF     = np.zeros((3,3,self.TipNode+1))   # axial reduction shape function (DOF by DOF by BldNodes+2)

        # integrate along the blade to find the generalized mass and stiffness, the twisted shape func and
        # the axial reduction shape function (blade z coordinate is coupled with x and y coordinates)

        # temporary variables
        TwstdSF             = np.zeros((2,3,2))
        TwstdSFOld          = np.zeros((2,3,2))
        AxRdBld             = np.zeros((3,3))
        AxRdBlforld         = np.zeros((3,3))

        for j in range(self.BldNodes):

            # generalized mass, ignore tip mass effect
            shapef1    = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldFl1Sh, 0)
            shapef2    = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldFl2Sh, 0)
            shapee     = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldEdgSh, 0)
            self.GenMassBF[0] = self.GenMassBF[0]+ self.BElmntMass[j]*shapef1*shapef1
            self.GenMassBF[1] = self.GenMassBF[1]+ self.BElmntMass[j]*shapef2*shapef2
            self.GenMassBE    = self.GenMassBE   + self.BElmntMass[j]*shapee * shapee

            # generalized stiffness
            shapef1 = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldFl1Sh, 2)
            shapef2 = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldFl2Sh, 2)
            shapee  = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldEdgSh, 2)
            self.GenStiffBF[0, 0] = self.GenStiffBF[0, 0] + self.StiffBF[j] * self.DRNodes[j]* shapef1 * shapef1
            self.GenStiffBF[0, 1] = self.GenStiffBF[0, 1] + self.StiffBF[j] * self.DRNodes[j]* shapef1 * shapef2
            self.GenStiffBF[1, 0] = self.GenStiffBF[1, 0] + self.StiffBF[j] * self.DRNodes[j]* shapef2 * shapef1
            self.GenStiffBF[1, 1] = self.GenStiffBF[1, 1] + self.StiffBF[j] * self.DRNodes[j]* shapef2 * shapef2
            self.GenStiffBE       = self.GenStiffBE+ self.StiffBE[j] * self.DRNodes[j]* shapee * shapee

            # twisted shape function
            self.TwistedSF[0, 0, j + 1, 2] = shapef1  * np.cos(self.Theta_s[j+1]*deg2rad)
            self.TwistedSF[1, 0, j + 1, 2] = -shapef1 * np.sin(self.Theta_s[j + 1]*deg2rad)
            self.TwistedSF[0, 1, j + 1, 2] = shapef2  * np.cos(self.Theta_s[j + 1]*deg2rad)
            self.TwistedSF[1, 1, j + 1, 2] = -shapef2 * np.sin(self.Theta_s[j + 1]*deg2rad)
            self.TwistedSF[0, 2, j + 1, 2] = shapee   * np.sin(self.Theta_s[j + 1]*deg2rad)
            self.TwistedSF[1, 2, j + 1, 2] = shapee   * np.cos(self.Theta_s[j + 1]*deg2rad)

            # the centrifugal stiffing term, this is actually not used in the ElastoDyn of OpenFAST 3.51
            #censtiff = FMomAbvNd[j]*self.DRNodes[j]*(TurbConfig['RotSpeed0']/60*2*np.pi)**2
            #shapef1  = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldFl1Sh, 1)
            #shapef2  = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldFl2Sh, 1)
            #shapee   = Shape(self.RNodesNorm[j], self.BldFlexL, self.BldEdgSh, 1)
            #self.CenStiffBF[0, 0] = self.CenStiffBF[0, 0] + censtiff * shapef1 * shapef1
            #self.CenStiffBF[1, 1] = self.CenStiffBF[1, 1] + censtiff * shapef2 * shapef2
            #self.CenStiffBE       = self.CenStiffBE + censtiff * shapee * shapee

            # Integrate to find the 1st derivative of the twisted shape functions:
            for I in range(2):
                for L in range(3):
                    TwstdSF[I, L, 1]               = self.TwistedSF[I, L, j + 1, 2]*0.5*self.DRNodes[j]
                    self.TwistedSF[I, L, j + 1, 1] = TwstdSF[I, L, 1]
            if j!=0:
                for I in range(2):
                    for L in range(3):
                        self.TwistedSF[I, L, j + 1, 1] = self.TwistedSF[I,L,j+1,1] + self.TwistedSF[I,L,j,1] + TwstdSFOld[I,L,1]

            # Integrate to find the zeroth derivative of the twisted shape functions:
            for I in range(2):
                for L in range(3):
                    TwstdSF[I, L, 0]               = self.TwistedSF[I, L, j + 1, 1]*0.5*self.DRNodes[j]
                    self.TwistedSF[I, L, j + 1, 0] = TwstdSF[I, L, 0]
            if j!=0:
                for I in range(2):
                    for L in range(3):
                        self.TwistedSF[I, L, j + 1, 0] = self.TwistedSF[I,L,j+1,0] + self.TwistedSF[I,L,j,0] + TwstdSFOld[I,L,0]

            # Integrate to find the blade axial reduction shape function
            for I in range(3):
                for L in range(3):
                    AxRdBld[I, L]                = 0.5*self.DRNodes[j]*(self.TwistedSF[0,I,j+1,1]*self.TwistedSF[0,L,j+1,1]+self.TwistedSF[1,I,j+1,1]*self.TwistedSF[1,L,j+1,1])
                    self.AxRedBldSF[I, L, j + 1] = AxRdBld[I, L]
            if j!=0:
                for I in range(3):
                    for L in range(3):
                        self.AxRedBldSF[I, L, j + 1] = self.AxRedBldSF[I,L,j+1]+self.AxRedBldSF[I,L,j]+AxRdBlforld[I,L]

            TwstdSFOld  = TwstdSF.copy()
            AxRdBlforld = AxRdBld.copy()

        # Calculate the mode natural frequencies which are needed for the generalized damping
        for I in range(2):
            self.FreqBF[I,0]  = np.sqrt(self.GenStiffBF[I,I]/( self.GenMassBF[I]))/np.pi/2
        self.FreqBE    = np.sqrt(self.GenStiffBE / (self.GenMassBE)) / np.pi/2

        BladeDamp = [BladeData['BldFlDmp1'],BladeData['BldFlDmp1']]
        for I in range(2):
            for L in range(2):
                self.GenDampBF[I,L] = 0.01*BladeDamp[I]*self.GenStiffBF[I,L]/(np.pi*self.FreqBF[I,0])
        self.GenDampBE              = 0.01 * BladeData['BldEdDmp1'] * self.GenStiffBE / (np.pi*self.FreqBE)

        # now, calculate the 2nd derivatives of the shape functions in lade roots and tip
        shapef1 = Shape(0.0, self.BldFlexL, self.BldFl1Sh, 2)
        shapef2 = Shape(0.0, self.BldFlexL, self.BldFl2Sh, 2)
        shapee  = Shape(0.0, self.BldFlexL, self.BldEdgSh, 2)
        self.TwistedSF[0, 0, 0, 2] = shapef1 * np.cos(self.Theta_s[0] * deg2rad)
        self.TwistedSF[1, 0, 0, 2] = -shapef1 * np.sin(self.Theta_s[0] * deg2rad)
        self.TwistedSF[0, 1, 0, 2] = shapef2 * np.cos(self.Theta_s[0] * deg2rad)
        self.TwistedSF[1, 1, 0, 2] = -shapef2 * np.sin(self.Theta_s[0] * deg2rad)
        self.TwistedSF[0, 2, 0, 2] = shapee * np.sin(self.Theta_s[0] * deg2rad)
        self.TwistedSF[1, 2, 0, 2] = shapee * np.cos(self.Theta_s[0] * deg2rad)

        shapef1 = Shape(1.0, self.BldFlexL, self.BldFl1Sh, 2)
        shapef2 = Shape(1.0, self.BldFlexL, self.BldFl2Sh, 2)
        shapee  = Shape(1.0, self.BldFlexL, self.BldEdgSh, 2)
        self.TwistedSF[0, 0, -1, 2] = shapef1 * np.cos(self.Theta_s[-1] * deg2rad)
        self.TwistedSF[1, 0, -1, 2] = -shapef1 * np.sin(self.Theta_s[-1] * deg2rad)
        self.TwistedSF[0, 1, -1, 2] = shapef2 * np.cos(self.Theta_s[-1] * deg2rad)
        self.TwistedSF[1, 1, -1, 2] = -shapef2 * np.sin(self.Theta_s[-1] * deg2rad)
        self.TwistedSF[0, 2, -1, 2] = shapee * np.sin(self.Theta_s[-1] * deg2rad)
        self.TwistedSF[1, 2, -1, 2] = shapee * np.cos(self.Theta_s[-1] * deg2rad)

        # integrate to get the first and zeroth order derivatives at the tip
        for I in range(2):
            for L in range(3):
                self.TwistedSF[I, L, -1, 1] = self.TwistedSF[I, L, self.BldNodes,1]  + TwstdSFOld[I,L,1]
                self.TwistedSF[I, L, -1, 0] = self.TwistedSF[I, L, self.BldNodes, 0] + TwstdSFOld[I, L, 0]

        self.TwistedSF[:,:, 0, 1]           = 0.0
        self.TwistedSF[:,:, 0, 0]           = 0.0
        self.AxRedBldSF[:,:,  0  ]          = 0.0
        for I in range(3):
            for L in range(3):
                self.AxRedBldSF[I,L, -1]    = self.AxRedBldSF[I,L, self.BldNodes]    + AxRdBlforld[I,L]

        self.AllBladeMass = self.NumBl*np.sum(self.BElmntMass)




class SetTowerPara:
    def __init__(self, TurbConfig,TowerData):
        self.TwrNodes  = TurbConfig['TwrNodes']
        self.TowerHt   = TurbConfig['TowerHt']
        self.TowerBsHt = TurbConfig['TowerBsHt']
        self.TwrFlexL  = TurbConfig['TowerHt'] - TurbConfig['TowerBsHt']

        self.DHNodes    = np.zeros((self.TwrNodes, 1))
        self.DHNodes[:] = self.TwrFlexL / self.TwrNodes
        self.HNodes     = np.zeros((self.TwrNodes, 1))
        self.HNodes[0]  = 0.5 * self.DHNodes[0]
        for J in range(1, self.TwrNodes): self.HNodes[J] = self.HNodes[J - 1] + 0.5 * ( self.DHNodes[J] + self.DHNodes[J - 1])
        self.HNodesNorm = self.HNodes / self.TwrFlexL
        self.TTopNode   = self.TwrNodes + 1

        # define tower cross-section properties
        self.MassT    = np.interp(self.HNodesNorm, TowerData['HtFract'], TowerData['TMassDen'])
        self.StiffTFA = np.interp(self.HNodesNorm, TowerData['HtFract'], TowerData['TwFAStif'])
        self.StiffTSS = np.interp(self.HNodesNorm, TowerData['HtFract'], TowerData['TwSSStif'])
        self.TwFAM1Sh = TowerData['TwFAM1Sh'] #shape functions
        self.TwFAM2Sh = TowerData['TwFAM2Sh'] #shape functions
        self.TwSSM1Sh = TowerData['TwSSM1Sh'] #shape functions
        self.TwSSM2Sh = TowerData['TwSSM2Sh'] #shape functions

        # tower top mass
        self.TwrTpMass = TurbConfig['AllBaldeMass'] + TurbConfig['HubMass'] + TurbConfig['NacMass'] + TurbConfig['YawBrMass']

        # tower element mass
        self.TElmntMass = np.zeros((self.TwrNodes, 1))
        for j in range(self.TwrNodes-1,-1,-1): self.TElmntMass[j] = self.MassT[j]*np.abs(self.DHNodes[j])

        # generalized mass and stiffness
        self.GenMassFA  = np.zeros((2, 2))  # generalized mass for tower fore-aft motion
        self.GenMassSS  = np.zeros((2, 2))  # generalized mass for blade side-side motion
        self.GenStiffFA = np.zeros((2, 2))  # generalized stiffness for tower fore-aft motion
        self.GenStiffSS = np.zeros((2, 2))  # generalized stiffness for tower side-side motion
        self.FreqTFA    = np.zeros((2, 1))  # natural frequency in fore-aft dir
        self.FreqTSS    = np.zeros((2, 1))  # natural frequency in side-side dir
        self.GenDampTFA  =  np.zeros((2, 2))  # generalized damping for tower fore-aft motion
        self.GenDampTSS  = np.zeros((2, 2))  # generalized damping for tower side-side motion

        for I in range(2): self.GenMassFA[I,I] = self.TwrTpMass; self.GenMassSS[I,I] = self.TwrTpMass

        # shape functions
        self.TwrFASF    = np.zeros((2, self.TTopNode+1, 3))  # DOF by DOF by derivative order
        self.TwrSSSF    = np.zeros((2, self.TTopNode+1, 3))
        self.AxRedTFA   = np.zeros((2, 2, self.TTopNode+1))
        self.AxRedTSS   = np.zeros((2, 2, self.TTopNode+1))

        # temporary variables
        AxRdTFA         = np.zeros((2, 2))
        AxRdTSS         = np.zeros((2, 2))
        AxRdTFAOld      = np.zeros((2, 2))
        AxRdTSSOld      = np.zeros((2, 2))

        for J in range(self.TwrNodes):
            self.TwrFASF[0, J + 1, 2] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwFAM1Sh, 2)
            self.TwrFASF[1, J + 1, 2] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwFAM2Sh, 2)
            self.TwrFASF[0, J + 1, 1] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwFAM1Sh, 1)
            self.TwrFASF[1, J + 1, 1] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwFAM2Sh, 1)
            self.TwrFASF[0, J + 1, 0] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwFAM1Sh, 0)
            self.TwrFASF[1, J + 1, 0] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwFAM2Sh, 0)

            self.TwrSSSF[0, J + 1, 2] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwSSM1Sh, 2)
            self.TwrSSSF[1, J + 1, 2] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwSSM2Sh, 2)
            self.TwrSSSF[0, J + 1, 1] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwSSM1Sh, 1)
            self.TwrSSSF[1, J + 1, 1] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwSSM2Sh, 1)
            self.TwrSSSF[0, J + 1, 0] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwSSM1Sh, 0)
            self.TwrSSSF[1, J + 1, 0] = Shape(self.HNodesNorm[J], self.TwrFlexL, self.TwSSM2Sh, 0)

            for I in range(2):
                self.GenMassFA[I,I] = self.GenMassFA[I,I]+self.TElmntMass[J]*self.TwrFASF[I,J+1,0]**2
                self.GenMassSS[I,I] = self.GenMassSS[I,I]+self.TElmntMass[J]*self.TwrSSSF[I,J+1,0]**2

            # integrate to find the generalized stiffness without gravity effect
            for I in range(2):
                for L in range(2):
                    self.GenStiffFA[I, L] += (self.StiffTFA[J]*self.DHNodes[J])* self.TwrFASF[I, J+1,2] * self.TwrFASF[L,J+1,2]
                    self.GenStiffSS[I, L] += (self.StiffTSS[J]*self.DHNodes[J])* self.TwrSSSF[I, J+1,2] * self.TwrSSSF[L,J+1,2]

            # integrate to find the axial shape reduction function
            for I in range(2):
                for L in range(2):
                    AxRdTFA[I, L] = 0.5 * self.DHNodes[J] * self.TwrFASF[I, J + 1, 1] * self.TwrFASF[L, J + 1, 1]
                    AxRdTSS[I, L] = 0.5 * self.DHNodes[J] * self.TwrSSSF[I, J + 1, 1] * self.TwrSSSF[L, J + 1, 1]
                    self.AxRedTFA[I, L, J + 1] = AxRdTFA[I, L]
                    self.AxRedTSS[I, L, J + 1] = AxRdTSS[I, L]
            if J!=0:
                for I in range(2):
                    for L in range(2):
                        self.AxRedTFA[I, L, J + 1] = self.AxRedTFA[I, L, J + 1] + self.AxRedTFA[I, L, J] + AxRdTFAOld[I,L]
                        self.AxRedTSS[I, L, J + 1] = self.AxRedTSS[I, L, J + 1] + self.AxRedTSS[I, L, J] + AxRdTSSOld[I,L]

            AxRdTFAOld = AxRdTFA.copy()
            AxRdTSSOld = AxRdTSS.copy()

        # calculate tower natural frequency without gravity effect
        for I in range(2):
            self.FreqTFA[I] = np.sqrt(self.GenStiffFA[I, I]/ (self.GenMassFA[I, I] - self.TwrTpMass))/2/np.pi
            self.FreqTSS[I] = np.sqrt(self.GenStiffSS[I, I]/ (self.GenMassSS[I, I] - self.TwrTpMass))/2/np.pi

        # calculate the generalized damping
        DampFA = [TowerData['TwrFADmp(1)'], TowerData['TwrFADmp(2)']]
        DampSS = [TowerData['TwrSSDmp(1)'], TowerData['TwrSSDmp(2)']]
        for I in range(2):
            for L in range(2):
                self.GenDampTFA[I, L] = (0.01 * DampFA[L]) * self.GenStiffFA[I, L] / (np.pi * self.FreqTFA[L])
                self.GenDampTSS[I, L] = (0.01 * DampSS[L]) * self.GenStiffSS[I, L] / (np.pi * self.FreqTSS[L])

        # tower shape functions at the top
        self.TwrFASF[0, self.TTopNode, 2] = Shape(1.0, self.TwrFlexL, self.TwFAM1Sh, 2)
        self.TwrFASF[1, self.TTopNode, 2] = Shape(1.0, self.TwrFlexL, self.TwFAM2Sh, 2)
        self.TwrFASF[0, self.TTopNode, 1] = Shape(1.0, self.TwrFlexL, self.TwFAM1Sh, 1)
        self.TwrFASF[1, self.TTopNode, 1] = Shape(1.0, self.TwrFlexL, self.TwFAM2Sh, 1)
        self.TwrFASF[0, self.TTopNode, 0] = Shape(1.0, self.TwrFlexL, self.TwFAM1Sh, 0)
        self.TwrFASF[1, self.TTopNode, 0] = Shape(1.0, self.TwrFlexL, self.TwFAM2Sh, 0)

        self.TwrSSSF[0, self.TTopNode, 2] = Shape(1.0, self.TwrFlexL, self.TwSSM1Sh, 2)
        self.TwrSSSF[1, self.TTopNode, 2] = Shape(1.0, self.TwrFlexL, self.TwSSM2Sh, 2)
        self.TwrSSSF[0, self.TTopNode, 1] = Shape(1.0, self.TwrFlexL, self.TwSSM1Sh, 1)
        self.TwrSSSF[1, self.TTopNode, 1] = Shape(1.0, self.TwrFlexL, self.TwSSM2Sh, 1)
        self.TwrSSSF[0, self.TTopNode, 0] = Shape(1.0, self.TwrFlexL, self.TwSSM1Sh, 0)
        self.TwrSSSF[1, self.TTopNode, 0] = Shape(1.0, self.TwrFlexL, self.TwSSM2Sh, 0)

        for I in range(2):
            for L in range(2):
                self.AxRedTFA[I,L,self.TTopNode] = self.AxRedTFA[I,L,self.TwrNodes]+ AxRdTFAOld[I,L]
                self.AxRedTSS[I,L,self.TTopNode] = self.AxRedTSS[I,L,self.TwrNodes]+ AxRdTSSOld[I,L]

        # tempdata = self.AxRedTSS.flatten('F')
        # for ii in range(tempdata.size):
        #     print(tempdata[ii])
        #
        # aa = 12


class SetDOFs:
    def __init__(self, DOF_Flag):
        self.DOF_Sg   =  1 # DOF index for platform surge
        self.DOF_Sw   =  2 # DOF index for platform sway
        self.DOF_Hv   =  3 # DOF index for platform heave
        self.DOF_R    =  4 # DOF index for platform roll
        self.DOF_P    =  5 # DOF index for platform pitch
        self.DOF_Y    =  6 # DOF index for platform yaw
        self.DOF_TFA1 =  7 # DOF index for 1st tower fore-aft mode
        self.DOF_TSS1 =  8 # DOF index for 1st tower side-to-side mode
        self.DOF_TFA2 =  9 # DOF index for 2nd tower fore-aft mode
        self.DOF_TSS2 = 10 # DOF index for 2nd tower side-to-side mode
        self.DOF_Yaw  = 11 # DOF index for nacelle-yaw
        self.DOF_RFrl = 12 # DOF index for rotor-furl, not actually used
        self.DOF_GeAz = 13 # DOF index for the generator azimuth
        self.DOF_DrTr = 14 # DOF index for drivetrain rotational-flexibility
        self.DOF_TFrl = 15 # DOF index for tail-furl, not actually used
        self.DOF_BF1  = [16, 19, 22]  # index for the first flap-wise mode for three blades
        self.DOF_BE   = [17, 20, 23]  # index for the first edge-wise mode for three blades
        self.DOF_BF2  = [18, 21, 24]  # index for the second flap-wise mode for three blades

        self.NactiveDOF  = 0
        self.ActiveDOFid = []
        self.ActiveDOFname = []
        # note that the rotor and tail furl DOFs are ignored
        self.FlexBlade = "false"
        self.FlexTower = "false"
        if DOF_Flag['PtfmSgDOF'].lower() == "true":
            self.NactiveDOF += 1
            self.ActiveDOFid.append(1)
            self.ActiveDOFname.append('Platform surge')
        if DOF_Flag['PtfmSwDOF'].lower()  == "true":
            self.NactiveDOF += 1
            self.ActiveDOFid.append(2)
            self.ActiveDOFname.append('Platform sway')
        if DOF_Flag['PtfmHvDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(3)
            self.ActiveDOFname.append('Platform heave')
        if DOF_Flag['PtfmRDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(4)
            self.ActiveDOFname.append('Platform roll, about x, Euler angles with sequence ZYX')
        if DOF_Flag['PtfmPDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(5)
            self.ActiveDOFname.append('Platform pitch, about y, Euler angles with sequence ZYX')
        if DOF_Flag['PtfmYDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(6)
            self.ActiveDOFname.append('Platform yaw, about z, Euler angles with sequence ZYX')
        if DOF_Flag['TwFADOF1'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(7)
            self.ActiveDOFname.append('First tower fore-aft mode')
            self.FlexTower = "true"
        if DOF_Flag['TwSSDOF1'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(8)
            self.ActiveDOFname.append('First tower side-side mode')
            self.FlexTower = "true"
        if DOF_Flag['TwFADOF2'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(9)
            self.ActiveDOFname.append('Second tower fore-aft mode')
            self.FlexTower = "true"
        if DOF_Flag['TwSSDOF2'].lower() == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(10)
            self.ActiveDOFname.append('Second tower side-side mode')
            self.FlexTower = "true"
        if DOF_Flag['YawDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(11)
            self.ActiveDOFname.append('Nacelle yaw mode')
        if DOF_Flag['GenDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(13)
            self.ActiveDOFname.append('Rotor azimuth mode')
        if DOF_Flag['DrTrDOF'].lower()  == "true" :
            self.NactiveDOF += 1
            self.ActiveDOFid.append(15)
            self.ActiveDOFname.append('Drivetrain torsional mode')
        if DOF_Flag['FlapDOF1'].lower()  == "true":
            self.NactiveDOF += 3
            self.ActiveDOFid.extend([16,19,22])
            self.ActiveDOFname.append('First blade flap-wise mode of blade 1')
            self.ActiveDOFname.append('First blade flap-wise mode of blade 2')
            self.ActiveDOFname.append('First blade flap-wise mode of blade 3')
            self.FlexBlade = "true"

        if DOF_Flag['FlapDOF2'].lower()  == "true":
            self.NactiveDOF += 3
            self.ActiveDOFid.extend([18, 21, 24])
            self.ActiveDOFname.append('Second blade flap-wise mode of blade 1')
            self.ActiveDOFname.append('Second blade flap-wise mode of blade 2')
            self.ActiveDOFname.append('Second blade flap-wise mode of blade 3')
            self.FlexBlade = "true"

        if DOF_Flag['EdgeDOF'].lower()  == "true":
            self.NactiveDOF += 3
            self.ActiveDOFid.extend([17, 20, 23])
            self.ActiveDOFname.append('First blade edge-wise mode of blade 1')
            self.ActiveDOFname.append('First blade edge-wise mode of blade 2')
            self.ActiveDOFname.append('First blade edge-wise mode of blade 3')
            self.FlexBlade = "true"


        # Zip and sort based on the DOF number
        self.ActiveDOFname = [x for _, x in sorted(zip(self.ActiveDOFid, self.ActiveDOFname))]
        self.ActiveDOFid.sort()



def SetMultibodyPara(TurbConfig,DOF_Flag,BladeData,TowerData):

    Blade      = SetBladePara(TurbConfig, BladeData)
    Platform   = SetPlatformPara(TurbConfig)
    TurbConfig['AllBaldeMass'] = Blade.AllBladeMass
    Tower      = SetTowerPara(TurbConfig,TowerData)
    Nacelle    = SetNacellePara(TurbConfig)
    Drivetrain = SetDriveTrainPara(TurbConfig)
    Hub        = SetHubPara(TurbConfig)
    DOFs       = SetDOFs(DOF_Flag)



    Parameters = SetParameters(Platform,Tower,Nacelle,Hub,Blade,Drivetrain,DOFs)
    return Parameters
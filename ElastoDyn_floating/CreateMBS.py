# create a multi-body system based on the...


import sympy as sm
import sympy.physics.mechanics as me
import numpy as np

deg2rad = 2 * np.pi / 360


class MBSsystem:
    def __init__(self, Parameters):

        # Specify points and reference frames from the bottom
        GlobalFrame = me.ReferenceFrame('G')  # inertial frame
        GlobalOrigin = me.Point('Go')

        PltFmFrame    = me.ReferenceFrame('P')
        PltFmPoint    = me.Point('Po')  # mass center point of platform
        PltFmRefPoint = me.Point('Pr') # reference point

        TowerFrame = []
        TowerPoint = []

        for i in range(Parameters.Tower.TwrNodes):
            TowerFrame.append(me.ReferenceFrame('T' + str(i + 1)))
            TowerPoint.append(me.Point('To' + str(i + 1)))

        TowerTopFrame = me.ReferenceFrame('TT')
        TowerTopPoint = me.Point('TTo')

        NacelleFrame = me.ReferenceFrame('N')
        NacellePoint = me.Point('No')

        TiltFrame = me.ReferenceFrame('TLT')  # tilted shaft coordinate system
        ApexPoint = me.Point('TLTo')  #

        HSSFrame = me.ReferenceFrame('HSS')  # high-speed shaft coordinate system
        HSSPoint = me.Point('HSSo')  #

        HubFrame = me.ReferenceFrame('H')  # also the low-speed shaft
        HubPoint = me.Point('Ho')

        BldConeFrame = []
        BldConePoint = []
        BldRootFrame = []
        BldRootPoint = []

        BladeFrame = [0] * int(Parameters.Blade.NumBl)
        BladePoint = [0] * int(Parameters.Blade.NumBl)

        for j in range(int(Parameters.Blade.NumBl)):
            BldConeFrame.append(me.ReferenceFrame('BC' + str(j + 1)))
            BldConePoint.append(me.Point('BC' + str(j + 1) + 'o'))
            BldRootFrame.append(me.ReferenceFrame('BR' + str(j + 1)))
            BldRootPoint.append(me.Point('BR' + str(j + 1) + 'o'))
            BladeFrame[j] = []
            BladePoint[j] = []

            for i in range(Parameters.Blade.BldNodes):
                BladeFrame[j].append(me.ReferenceFrame('B' + str(j + 1) + 'E' + str(i + 1)))
                BladePoint[j].append(me.Point('B' + str(j + 1) + 'Eo' + str(i + 1)))

        # Configurate DOFs
        GenCoord = []
        GenSpeed = []
        GenAcc = []
        ACCtoReplace = []
        StoReplace = []

        t = me.dynamicsymbols._t

        ReplaceString = ''
        for i in range(24):
            if i + 1 in Parameters.DOFs.ActiveDOFid:
                GenCoord.append(me.dynamicsymbols('q' + str(i + 1)))
                # GenSpeed.append(me.dynamicsymbols('u' + str(i + 1)))
                GenSpeed.append(GenCoord[i].diff(t))
                GenAcc.append(me.dynamicsymbols('a' + str(i + 1)))
                StoReplace.append(GenCoord[i].diff(t, 1))
                ACCtoReplace.append(GenCoord[i].diff(t, 2))
                # tempString = 'GenCoord[' + str(i) + '].diff(t,2): GenAcc[' + str(i) + '],'
                # ReplaceString = ReplaceString + tempString
            else:
                exec('GenCoord.append(0)')
                exec('GenSpeed.append(0)')
                exec('GenAcc.append(0)')

        # qdot2repl = {' + ReplaceString + '}

        ActiveDOF = [x for x in GenCoord if x != 0]
        ActiveGenSpeed = [x for x in GenSpeed if x != 0]

        self.ActiveDOF = ActiveDOF
        self.ActiveGenSpeed = ActiveGenSpeed
        self.ActiveAcc = [x for x in GenAcc if x != 0]

        qdot2_repl = dict(zip(ACCtoReplace, self.ActiveAcc))
        qdot_repl = dict(zip(StoReplace, self.ActiveGenSpeed))

        # def Numerify(Partial):
        #     func = sm.Matrix([0, 0, 0])
        #     func[0] = Partial.dot(GlobalFrame.x)
        #     func[1] = Partial.dot(GlobalFrame.y)
        #     func[2] = Partial.dot(GlobalFrame.z)
        #
        #     return sm.lambdify((ActiveDOF, ActiveGenSpeed), [func])

        # get some parameters
        DOFS = Parameters.DOFs
        TwrFASF = Parameters.Tower.TwrFASF
        TwrSSSF = Parameters.Tower.TwrSSSF
        AxRedTFA = Parameters.Tower.AxRedTFA
        AxRedTSS = Parameters.Tower.AxRedTSS
        Hnodes = Parameters.Tower.HNodes
        TTopNode = int(Parameters.Tower.TTopNode)

        PtfmCM_xyz = Parameters.Platform.PtfmCM_xyz

        # set rotational and translational relationship between platform and inertial frames
        PltFmFrame.orient_body_fixed(GlobalFrame,
                                     (GenCoord[DOFS.DOF_Y - 1], GenCoord[DOFS.DOF_P - 1], GenCoord[DOFS.DOF_R - 1]),
                                     'zyx')

        PlatFormAngvec = 0
        if GenCoord[DOFS.DOF_Y - 1] !=0: PlatFormAngvec += GenCoord[DOFS.DOF_Y - 1].diff(t) * GlobalFrame.z
        if GenCoord[DOFS.DOF_P - 1] !=0: PlatFormAngvec += GenCoord[DOFS.DOF_P - 1].diff(t) * GlobalFrame.y
        if GenCoord[DOFS.DOF_R - 1] !=0: PlatFormAngvec += GenCoord[DOFS.DOF_R - 1].diff(t) * GlobalFrame.x

        PltFmFrame.set_ang_vel(GlobalFrame,PlatFormAngvec)

        PltFmRefPoint.set_pos(GlobalOrigin,
                           GenCoord[DOFS.DOF_Sg - 1]  * GlobalFrame.x + GenCoord[DOFS.DOF_Sw - 1]  * GlobalFrame.y +
                           GenCoord[DOFS.DOF_Hv - 1]  * GlobalFrame.z)
        PltFmPoint.set_pos(PltFmRefPoint,
                           PtfmCM_xyz[0] * PltFmFrame.x + PtfmCM_xyz[1] * PltFmFrame.y +
                            + PtfmCM_xyz[2] * PltFmFrame.z)

        GlobalOrigin.set_vel(GlobalFrame, 0)
        PltFmPoint.set_vel(PltFmFrame, 0)
        PltFmRefPoint.set_vel(PltFmFrame, 0)
        PltFmPoint.v2pt_theory(PltFmRefPoint, GlobalFrame, PltFmFrame)

        self.PltFm_acc = PltFmPoint.acc(GlobalFrame).xreplace(qdot2_repl)
        self.PltFm_ang_acc = PltFmFrame.ang_acc_in(GlobalFrame).xreplace(qdot2_repl)

        PltFm_P_vel     = []
        PltFmRef_P_vel  = []
        PltFm_P_ang_vel = []

        for dq in ActiveGenSpeed:
            PltFm_P_vel.append(PltFmPoint.vel(GlobalFrame).diff(dq, GlobalFrame).express(GlobalFrame))
            PltFmRef_P_vel.append(PltFmRefPoint.vel(GlobalFrame).diff(dq, GlobalFrame).express(GlobalFrame))
            PltFm_P_ang_vel.append(PltFmFrame.ang_vel_in(GlobalFrame).diff(dq, GlobalFrame).express(GlobalFrame))

        self.PltFm_P_vel = PltFm_P_vel
        self.PltFmRef_P_vel = PltFmRef_P_vel
        self.PltFm_P_ang_vel = PltFm_P_ang_vel

        # set rotational and translational relationship between tower and platform

        TowerEle_acc = []

        TowerEle_ang_acc = []
        TowerEle_P_vel = [0] * Parameters.Tower.TwrNodes
        TowerEle_P_ang_vel = [0] * Parameters.Tower.TwrNodes

        for i in range(Parameters.Tower.TwrNodes):  # tower element frames
            ThetaFA = TwrFASF[0, i + 1, 1] * GenCoord[DOFS.DOF_TFA1 - 1] + TwrFASF[1, i + 1, 1] * GenCoord[
                DOFS.DOF_TFA2 - 1]
            ThetaSS = -TwrSSSF[0, i + 1, 1] * GenCoord[DOFS.DOF_TSS1 - 1] - TwrSSSF[1, i + 1, 1] * GenCoord[
                DOFS.DOF_TSS2 - 1]
            TowerFrame[i].orient_body_fixed(PltFmFrame, (0, ThetaFA, ThetaSS), 'zyx')

            if ThetaFA == 0 and ThetaSS != 0:
                TowerFrame[i].set_ang_vel(PltFmFrame, ThetaSS.diff(t) * PltFmFrame.x)
            elif ThetaSS == 0 and ThetaFA != 0:
                TowerFrame[i].set_ang_vel(PltFmFrame, ThetaFA.diff(t) * PltFmFrame.y)
            elif ThetaSS != 0 and ThetaFA != 0:
                TowerFrame[i].set_ang_vel(PltFmFrame, ThetaFA.diff(t) * PltFmFrame.y + ThetaSS.diff(t) * PltFmFrame.x)

            # TowerFrame[i].set_ang_vel(GlobalFrame, TowerFrame[i].ang_vel_in(GlobalFrame).xreplace(qdot_repl))

            Delta_x = TwrFASF[0, i + 1, 0] * GenCoord[DOFS.DOF_TFA1 - 1] + TwrFASF[1, i + 1, 0] * GenCoord[
                DOFS.DOF_TFA2 - 1]
            Delta_y = TwrSSSF[0, i + 1, 0] * GenCoord[DOFS.DOF_TSS1 - 1] + TwrSSSF[1, i + 1, 0] * GenCoord[
                DOFS.DOF_TSS2 - 1]
            Delta_z = float(Parameters.Tower.HNodes[i]+Parameters.Tower.TowerBsHt-PtfmCM_xyz[2]) - 0.5 * (
                    AxRedTFA[0, 0, i + 1] * GenCoord[DOFS.DOF_TFA1 - 1] ** 2
                    + AxRedTFA[1, 1, i + 1] * GenCoord[DOFS.DOF_TFA2 - 1] ** 2
                    + 2.0 * AxRedTFA[0, 1, i + 1] * GenCoord[DOFS.DOF_TFA1 - 1] * GenCoord[DOFS.DOF_TFA2 - 1]
                    + AxRedTSS[0, 0, i + 1] * GenCoord[DOFS.DOF_TSS1 - 1] ** 2
                    + AxRedTSS[1, 1, i + 1] * GenCoord[DOFS.DOF_TSS2 - 1] ** 2
                    + 2.0 * AxRedTSS[0, 1, i + 1] * GenCoord[DOFS.DOF_TSS1 - 1] * GenCoord[DOFS.DOF_TSS2 - 1])

            TowerPoint[i].set_pos(PltFmPoint, Delta_x * PltFmFrame.x + Delta_y * PltFmFrame.y + Delta_z * PltFmFrame.z)
            TowerPoint[i].v1pt_theory(PltFmPoint, GlobalFrame, PltFmFrame)
            TowerEle_acc.append(TowerPoint[i].acc(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl))
            TowerEle_ang_acc.append(TowerFrame[i].ang_acc_in(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl))

            TowerEle_P_vel[i] = []
            TowerEle_P_ang_vel[i] = []

            for q in ActiveGenSpeed:
                TowerEle_P_vel[i].append(TowerPoint[i].vel(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))
                TowerEle_P_ang_vel[i].append(TowerFrame[i].ang_vel_in(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))

        self.TowerEle_acc = TowerEle_acc
        self.TowerEle_ang_acc = TowerEle_ang_acc
        self.TowerEle_P_vel = TowerEle_P_vel
        self.TowerEle_P_ang_vel = TowerEle_P_ang_vel

        ThetaFA = TwrFASF[0, -1, 1] * GenCoord[DOFS.DOF_TFA1 - 1] + TwrFASF[1, -1, 1] * GenCoord[DOFS.DOF_TFA2 - 1]
        ThetaSS = -TwrSSSF[0, -1, 1] * GenCoord[DOFS.DOF_TSS1 - 1] -TwrSSSF[1, -1, 1] * GenCoord[DOFS.DOF_TSS2 - 1]
        # ^ right hand rule coord: positive displace in y is negative rotation around x

        TowerTopFrame.orient_body_fixed(PltFmFrame, (0, ThetaFA, ThetaSS), 'zyx')  # tower top frame
        if ThetaFA == 0 and ThetaSS != 0:
            TowerTopFrame.set_ang_vel(PltFmFrame, ThetaSS.diff(t) * PltFmFrame.x)
        elif ThetaSS == 0 and ThetaFA != 0:
            TowerTopFrame.set_ang_vel(PltFmFrame, ThetaFA.diff(t) * PltFmFrame.y)
        elif ThetaSS != 0 and ThetaFA != 0:
            TowerTopFrame.set_ang_vel(PltFmFrame, ThetaFA.diff(t) * PltFmFrame.y + ThetaSS.diff(t) * PltFmFrame.x)

        # TowerTopFrame.set_ang_vel(GlobalFrame, TowerTopFrame.ang_vel_in(GlobalFrame).xreplace(qdot_repl))

        Delta_x = 1.0 * GenCoord[DOFS.DOF_TFA1 - 1] + 1.0 * GenCoord[DOFS.DOF_TFA2 - 1]
        Delta_y = 1.0 * GenCoord[DOFS.DOF_TSS1 - 1] + 1.0 * GenCoord[DOFS.DOF_TSS2 - 1]
        Delta_z = float(Parameters.Tower.TowerHt -PtfmCM_xyz[2]) - 0.5 * (
                AxRedTFA[0, 0, TTopNode] * GenCoord[DOFS.DOF_TFA1 - 1] ** 2  # check tip node later
                + AxRedTFA[1, 1, TTopNode] * GenCoord[DOFS.DOF_TFA2 - 1] ** 2
                + 2.0 * AxRedTFA[0, 1, TTopNode] * GenCoord[DOFS.DOF_TFA1 - 1] * GenCoord[DOFS.DOF_TFA2 - 1]
                + AxRedTSS[0, 0, TTopNode] * GenCoord[DOFS.DOF_TSS1 - 1] ** 2
                + AxRedTSS[1, 1, TTopNode] * GenCoord[DOFS.DOF_TSS2 - 1] ** 2
                + 2.0 * AxRedTSS[0, 1, TTopNode] * GenCoord[DOFS.DOF_TSS1 - 1] * GenCoord[DOFS.DOF_TSS2 - 1])
        TowerTopPoint.set_pos(PltFmPoint, Delta_x * PltFmFrame.x + Delta_y * PltFmFrame.y + Delta_z * PltFmFrame.z)
        TowerTopPoint.set_vel(TowerTopFrame, 0)
        TowerTopPoint.v1pt_theory(PltFmPoint, GlobalFrame, PltFmFrame)

        self.TowerTop_acc = TowerTopPoint.acc(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)
        self.TowerTop_ang_acc = TowerTopFrame.ang_acc_in(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)

        TowerTop_P_vel = []
        TowerTop_P_ang_vel = []

        for q in ActiveGenSpeed:
            TowerTop_P_vel.append(TowerTopPoint.vel(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))
            TowerTop_P_ang_vel.append(TowerTopFrame.ang_vel_in(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))

        self.TowerTop_P_vel = TowerTop_P_vel
        self.TowerTop_P_ang_vel = TowerTop_P_ang_vel

        # set rotational and translational relationships for nacelle, shaft, hub and blade roots

        NacelleFrame.orient_axis(TowerTopFrame, GenCoord[DOFS.DOF_Yaw - 1], TowerTopFrame.z)  # nacelle frame

        Delta_x = float(Parameters.Nacelle.NacCM_xyz[0])
        Delta_y = float(Parameters.Nacelle.NacCM_xyz[1])
        Delta_z = float(Parameters.Nacelle.NacCM_xyz[2])
        NacellePoint.set_pos(TowerTopPoint,
                             Delta_x * NacelleFrame.x + Delta_y * NacelleFrame.y + Delta_z * TowerTopFrame.z)
        NacellePoint.v1pt_theory(TowerTopPoint, GlobalFrame, TowerTopFrame)

        self.Nacelle_acc = NacellePoint.acc(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)
        self.Nacelle_ang_acc = NacelleFrame.ang_acc_in(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)

        Nacelle_P_vel = []
        Nacelle_P_ang_vel = []

        for q in ActiveGenSpeed:
            Nacelle_P_vel.append(NacellePoint.vel(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))
            Nacelle_P_ang_vel.append(NacelleFrame.ang_vel_in(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))

        self.Nacelle_P_vel = Nacelle_P_vel
        self.Nacelle_P_ang_vel = Nacelle_P_ang_vel
        # In[27]:

        TiltFrame.orient_axis(NacelleFrame, -Parameters.DriveTrain.ShftTilt * deg2rad, NacelleFrame.y)
        # tilted shaft, note that the angle is defined negative if rotor tilts up

        Twr2Shft = float(Parameters.DriveTrain.Twr2Shft)
        OverHang = float(Parameters.Hub.OverHang)
        HubCm = float(Parameters.Hub.HubCM)

        ApexPoint.set_pos(TowerTopPoint, Twr2Shft * NacelleFrame.z + OverHang * TiltFrame.x)
        ApexPoint.set_vel(NacelleFrame, 0)

        HSSFrame.orient_axis(TiltFrame, GenCoord[DOFS.DOF_GeAz - 1], TiltFrame.x)

        HSSPoint.set_pos(TowerTopPoint, Twr2Shft * NacelleFrame.z)  # the HSS is not implemented in openfast
        HSSPoint.set_vel(HSSFrame, 0)
        HSSPoint.v1pt_theory(TowerTopPoint, GlobalFrame, TowerTopFrame).xreplace(qdot2_repl)

        self.HSS_acc = HSSPoint.acc(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)
        self.HSS_ang_acc = HSSFrame.ang_acc_in(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)

        HSS_P_vel = []
        HSS_P_ang_vel = []

        for q in ActiveGenSpeed:
            HSS_P_vel.append(HSSPoint.vel(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))
            HSS_P_ang_vel.append(HSSFrame.ang_vel_in(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))

        self.HSS_P_vel = HSS_P_vel
        self.HSS_P_ang_vel = HSS_P_ang_vel

        HubFrame.orient_axis(HSSFrame, GenCoord[DOFS.DOF_DrTr - 1], HSSFrame.x)

        HubPoint.set_pos(ApexPoint, HubCm * TiltFrame.x)
        HubPoint.set_vel(NacelleFrame, 0)
        HubPoint.v1pt_theory(TowerTopPoint, GlobalFrame, TowerTopFrame)

        self.Hub_acc = HubPoint.acc(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)
        self.Hub_ang_acc = HubFrame.ang_acc_in(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)

        Hub_P_vel = []
        Hub_P_ang_vel = []

        for q in ActiveGenSpeed:
            Hub_P_vel.append(HubPoint.vel(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))
            Hub_P_ang_vel.append(HubFrame.ang_vel_in(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))

        self.Hub_P_vel = Hub_P_vel
        self.Hub_P_ang_vel = Hub_P_ang_vel

        HubRadius = float(Parameters.Hub.HubRad)

        ## assume zero pitch for now, they should depend on the pitch actuator dynamics, and adds three more
        # DOFs, in Elastodyn, these dynamics are not modelled
        beta = [0, 0, 0]

        for j in range(int(Parameters.Blade.NumBl)):
            BldConeFrame[j].orient_body_fixed(HubFrame, (j * 2 / 3 * np.pi, Parameters.Blade.PreCone * deg2rad, 0),
                                              'xyz')  # blade1 cone
            BldConePoint[j].set_pos(ApexPoint, HubRadius * BldConeFrame[j].z)
            BldConePoint[j].set_vel(HubFrame, 0)
            BldConePoint[j].v1pt_theory(HubPoint, GlobalFrame, HubFrame)
            BldRootFrame[j].orient_axis(BldConeFrame[j], beta[j], BldConeFrame[j].z)
            BldRootPoint[j].set_pos(BldConePoint[j], 0)
            BldRootPoint[j].set_vel(HubFrame, 0)  # assume the root is stiffly connected to hub
            BldRootPoint[j].v1pt_theory(HubPoint, GlobalFrame, HubFrame)

        # blade elemments
        TwistedSF = Parameters.Blade.TwistedSF
        RNodes = Parameters.Blade.RNodes
        AxRedBldSF = Parameters.Blade.AxRedBldSF
        for i in range(Parameters.Blade.BldNodes):  # Parameters.Blade.BldNodes):
            for j in range(int(Parameters.Blade.NumBl)):
                Theta_S = float(Parameters.Blade.Theta_s[i] * deg2rad)
                ThetaOoP = TwistedSF[0, 0, i + 1, 1] * GenCoord[DOFS.DOF_BF1[j] - 1] \
                           + TwistedSF[0, 1, i + 1, 1] * GenCoord[DOFS.DOF_BF2[j] - 1] \
                           + TwistedSF[0, 2, i + 1, 1] * GenCoord[DOFS.DOF_BE[j] - 1]
                ThetaIP = -TwistedSF[1, 0, i + 1, 1] * GenCoord[DOFS.DOF_BF1[j] - 1] \
                          - TwistedSF[1, 1, i + 1, 1] * GenCoord[DOFS.DOF_BF2[j] - 1] \
                          - TwistedSF[1, 2, i + 1, 1] * GenCoord[DOFS.DOF_BE[j] - 1]
                ThetaLxb = float(np.cos(Theta_S)) * ThetaIP - float(np.sin(Theta_S)) * ThetaOoP
                ThetaLyb = float(np.sin(Theta_S)) * ThetaIP + float(np.cos(Theta_S)) * ThetaOoP

                BladeFrame[j][i].orient_body_fixed(BldRootFrame[j], (0, ThetaOoP, ThetaIP), 'zyx')
                if ThetaIP == 0 and ThetaOoP != 0:
                    BladeFrame[j][i].set_ang_vel(BldRootFrame[j], ThetaOoP.diff(t) * BldRootFrame[j].y)
                elif ThetaOoP == 0 and ThetaIP != 0:
                    BladeFrame[j][i].set_ang_vel(BldRootFrame[j], ThetaIP.diff(t) * BldRootFrame[j].x)
                elif ThetaOoP != 0 and ThetaIP != 0:
                    BladeFrame[j][i].set_ang_vel(BldRootFrame[j],
                                                 ThetaIP.diff(t) * BldRootFrame[j].x + ThetaOoP.diff(t) * BldRootFrame[
                                                     j].y)

                Delta_x = TwistedSF[0, 0, i + 1, 0] * GenCoord[DOFS.DOF_BF1[j] - 1] \
                          + TwistedSF[0, 1, i + 1, 0] * GenCoord[DOFS.DOF_BF2[j] - 1] \
                          + TwistedSF[0, 2, i + 1, 0] * GenCoord[DOFS.DOF_BE[j] - 1]

                Delta_y = TwistedSF[1, 0, i + 1, 0] * GenCoord[DOFS.DOF_BF1[j] - 1] \
                          + TwistedSF[1, 1, i + 1, 0] * GenCoord[DOFS.DOF_BF2[j] - 1] \
                          + TwistedSF[1, 2, i + 1, 0] * GenCoord[DOFS.DOF_BE[j] - 1]

                Delta_z = float(RNodes[i]) - 0.5 * (AxRedBldSF[0, 0, i + 1] * GenCoord[DOFS.DOF_BF1[j] - 1] ** 2 \
                                                    + AxRedBldSF[1, 1, i + 1] * GenCoord[DOFS.DOF_BF2[j] - 1] ** 2 \
                                                    + AxRedBldSF[2, 2, i + 1] * GenCoord[DOFS.DOF_BE[j] - 1] ** 2 \
                                                    + 2.0 * AxRedBldSF[0, 1, i + 1] * GenCoord[DOFS.DOF_BF1[j] - 1] *
                                                    GenCoord[DOFS.DOF_BF2[j] - 1] \
                                                    + 2.0 * AxRedBldSF[1, 2, i + 1] * GenCoord[DOFS.DOF_BF2[j] - 1] *
                                                    GenCoord[DOFS.DOF_BE[j] - 1] \
                                                    + 2.0 * AxRedBldSF[0, 2, i + 1] * GenCoord[DOFS.DOF_BF1[j] - 1] *
                                                    GenCoord[DOFS.DOF_BE[j] - 1])

                BladePoint[j][i].set_pos(BldRootPoint[j],
                                         Delta_x * BldRootFrame[j].x + Delta_y * BldRootFrame[j].y + Delta_z *
                                         BldRootFrame[j].z)
                BladePoint[j][i].v1pt_theory(HubPoint, GlobalFrame, HubFrame)

        BladeEle_acc = [[0 for col in range(Parameters.Blade.BldNodes)] for col in range(int(Parameters.Blade.NumBl))]
        BladeEle_ang_acc = [[0 for col in range(Parameters.Blade.BldNodes)] for col in
                            range(int(Parameters.Blade.NumBl))]
        BladeEle_P_vel = [[0 for col in range(Parameters.Blade.BldNodes)] for col in range(int(Parameters.Blade.NumBl))]
        BladeEle_P_ang_vel = [[0 for col in range(Parameters.Blade.BldNodes)] for col in
                              range(int(Parameters.Blade.NumBl))]

        for i in range(Parameters.Blade.BldNodes):  # Parameters.Blade.BldNodes):
            for j in range(int(Parameters.Blade.NumBl)):
                BladeEle_acc[j][i] = BladePoint[j][i].acc(GlobalFrame).express(GlobalFrame).xreplace(qdot2_repl)
                BladeEle_ang_acc[j][i] = BladeFrame[j][i].ang_acc_in(GlobalFrame).express(GlobalFrame).xreplace(
                    qdot2_repl)

                BladeEle_P_vel[j][i] = []
                BladeEle_P_ang_vel[j][i] = []

                for q in ActiveGenSpeed:
                    BladeEle_P_vel[j][i].append(
                        BladePoint[j][i].vel(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))
                    BladeEle_P_ang_vel[j][i].append(
                        BladeFrame[j][i].ang_vel_in(GlobalFrame).diff(q, GlobalFrame).express(GlobalFrame))

        self.BladeEle_acc = BladeEle_acc
        self.BladeEle_ang_acc = BladeEle_ang_acc
        self.BladeEle_P_vel = BladeEle_P_vel
        self.BladeEle_P_ang_vel = BladeEle_P_ang_vel

        self.GlobalFrame = GlobalFrame
        self.HubFrame = HubFrame
        self.NacelleFrame = NacelleFrame
        self.HSSFrame = HSSFrame
        self.TowerTopFrame = TowerTopFrame
        self.PltFmFrame = PltFmFrame
        self.DCM_TLT2G  = GlobalFrame.dcm(TiltFrame)

    def NumerifyMBSv2(self, Parameters):

        # convert the symbolic MBS to numerical
        GlobalFrame = self.GlobalFrame
        ActiveDOF = self.ActiveDOF
        ActiveGenSpeed = self.ActiveGenSpeed
        ActiveAcc = self.ActiveAcc

        # firstly, Calculate the generalized mass matrix----------------
        Frs_bar = []
        dM = []
        for i in range(len(ActiveDOF)):
            dM.append([0] * len(ActiveDOF))

        for k in range(len(ActiveDOF)):

            Frs = 0

            # firstly the blade contribution--------------------
            for i in range(Parameters.Blade.BldNodes):  # Parameters.Blade.BldNodes
                masses = float(Parameters.Blade.BElmntMass[i])

                for j in range(int(Parameters.Blade.NumBl)):
                    vr = self.BladeEle_P_vel[j][i][k]
                    Rs = -masses * self.BladeEle_acc[j][i]
                    Frs += vr.dot(Rs)

                    # in Elastodyn the central inertia dyadic of blade element is not defined, thus they are ignored
                    # wr = self.BladeEle_P_ang_vel[j][i][k]
                    # Ts = -(self.BladeEle_ang_acc[j][i].dot(inertias) +
                    #        me.cross(self.BladeEle_ang_acc[j][i], inertias).dot(self.BladeEle_ang_acc[j][i]))
                    # Frs += wr.dot(Ts)
            # Jacobian
            for m in range(k + 1):
                dM[k][m] = Frs.diff(ActiveAcc[m])

            # Hub and low speed shaft----------------------
            masses = float(Parameters.Hub.HubMass)  # the lss shaft mass is included in the nacelle
            inertias = float(Parameters.Hub.HubIner) * me.outer(self.HubFrame.x,
                                                                self.HubFrame.x)  # float(Parameters.Hub.HubIner)*HubFrame.x

            vr = self.Hub_P_vel[k]
            Rs = -masses * self.Hub_acc
            Fr_temp = vr.dot(Rs)
            Frs += Fr_temp

            # Jacobian
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            wr = self.Hub_P_ang_vel[k]
            Ts = -(self.Hub_ang_acc.dot(inertias) +
                   me.cross(self.Hub_ang_acc, inertias).dot(self.Hub_ang_acc))
            Fr_temp = wr.dot(Ts)
            Frs += Fr_temp
            # Jacobian
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            # High speed shaft-------------------------
            # masses   = float(Parameters.Hub.HubMass)  # the hss shaft mass is included in the nacelle
            inertias = float(Parameters.DriveTrain.GeneIner) * me.outer(self.HSSFrame.x,
                                                                        self.HSSFrame.x)  # float(Parameters.Hub.HubIner)*HubFrame.x

            wr = self.HSS_P_ang_vel[k]
            Ts = -(self.HSS_ang_acc.dot(inertias) +
                   me.cross(self.HSS_ang_acc, inertias).dot(self.HSS_ang_acc))
            Fr_temp = wr.dot(Ts)
            Frs += Fr_temp
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])


            # Nacelle-----------------------------------
            masses = float(Parameters.Nacelle.NacMass)  #
            inertias = float(Parameters.Nacelle.NacYIner) * me.outer(self.NacelleFrame.z, self.NacelleFrame.z)
            # the nacelle central inertia dyadic in x and y are not defined

            vr = self.Nacelle_P_vel[k]
            Rs = -masses * self.Nacelle_acc
            Fr_temp = vr.dot(Rs)
            Frs += Fr_temp
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])


            wr = self.Nacelle_P_ang_vel[k]
            Ts = -(self.Nacelle_ang_acc.dot(inertias) +
                   me.cross(self.Nacelle_ang_acc, inertias).dot(self.Nacelle_ang_acc))
            Fr_temp = wr.dot(Ts)
            Frs += Fr_temp
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            # Tower top----------------------------------
            masses = float(Parameters.Nacelle.YawBrMass)

            vr = self.TowerTop_P_vel[k]
            Rs = -masses * self.TowerTop_acc
            Fr_temp = vr.dot(Rs)
            Frs += Fr_temp
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            # Tower elements----------------------------
            for i in range(int(Parameters.Tower.TwrNodes)):
                masses = float(Parameters.Tower.TElmntMass[i])
                vr = self.TowerEle_P_vel[i][k]
                Rs = -masses * self.TowerEle_acc[i]
                Fr_temp = vr.dot(Rs)
                Frs += Fr_temp
                for m in range(k + 1):
                    dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            # Platform --------------------------------
            masses = float(Parameters.Platform.PtfmMass)  #
            inertias = float(Parameters.Platform.PtfmRIner) * me.outer(self.PltFmFrame.x, self.PltFmFrame.x) \
                       + float(Parameters.Platform.PtfmPIner) * me.outer(self.PltFmFrame.y, self.PltFmFrame.y) \
                       + float(Parameters.Platform.PtfmYIner) * me.outer(self.PltFmFrame.z, self.PltFmFrame.z)

            vr = self.PltFm_P_vel[k]
            Rs = -masses * self.PltFm_acc
            Fr_temp = vr.dot(Rs)
            Frs += Fr_temp
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            wr = self.PltFm_P_ang_vel[k]
            Ts = -(self.PltFm_ang_acc.dot(inertias) +
                   me.cross(self.PltFm_ang_acc, inertias).dot(self.PltFm_ang_acc))
            Fr_temp = wr.dot(Ts)
            Frs += Fr_temp
            for m in range(k + 1):
                dM[k][m] += Fr_temp.diff(ActiveAcc[m])

            # Fr_bar.append(Fr)
            Frs_bar.append(Frs)

        # calculate the Taylor expansion of Frs
        Frs = sm.Matrix(Frs_bar)
        a_zerod = {x: 0 for x in ActiveAcc}

        # M_de = Frs.jacobian(ActiveAcc).xreplace(a_zerod)  # in principle, the acc terms are linear summation
        # the jacobian should give a matrix free of acc terms
        # while, in current sympy. Some extremely small values present
        # as mulplification factors before part of the acc terms,
        # which can be a numerically bug in current sympy.
        # xreplace is used here to avoid this issue.

        for k in range(len(ActiveDOF)):
            for m in range(k+1,len(ActiveDOF)):
                dM[k][m] = dM[m][k]
        M_de   = sm.Matrix(dM)
        self.M = sm.lambdify((ActiveDOF, ActiveGenSpeed), M_de.xreplace(a_zerod))

        #q_vals = np.array([0, 0])
        #qdt_vals = np.array([0,0.598996996879578])

        self.M0 = sm.lambdify((ActiveDOF, ActiveGenSpeed), Frs.xreplace(a_zerod))
        self.Frs_bar = Frs_bar
        self.M_de = M_de
        self.dM   = dM


        DCM = [[0 for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                DCM[i][j]=self.DCM_TLT2G[i+j]

        self.DCM_TLT2G_num = sm.lambdify((ActiveDOF, ActiveGenSpeed), sm.Matrix(DCM), 'numpy') # The generalized sppeds are dummy values


        # numerify all the partial velocity and partial angular velocity

        def Numerify(Partial):
            func = sm.Matrix([0, 0, 0])
            func[0] = Partial.dot(GlobalFrame.x)
            func[1] = Partial.dot(GlobalFrame.y)
            func[2] = Partial.dot(GlobalFrame.z)

            return sm.lambdify((ActiveDOF, ActiveGenSpeed), [func], 'numpy')

        self.BladeEle_P_vel_num = create_same_shape_list(self.BladeEle_P_vel)
        self.BladeEle_P_ang_vel_num = create_same_shape_list(self.BladeEle_P_ang_vel)

        for i in range(len(self.BladeEle_P_vel_num)):
            for j in range(len(self.BladeEle_P_vel_num[0])):
                for k in range(len(self.BladeEle_P_vel_num[0][0])):
                    self.BladeEle_P_vel_num[i][j][k] = Numerify(self.BladeEle_P_vel[i][j][k])
                    self.BladeEle_P_ang_vel_num[i][j][k] = Numerify(self.BladeEle_P_ang_vel[i][j][k])

        self.Hub_P_vel_num = create_same_shape_list(self.Hub_P_vel)
        self.Hub_P_ang_vel_num = create_same_shape_list(self.Hub_P_ang_vel)
        self.Nacelle_P_vel_num = create_same_shape_list(self.Nacelle_P_vel)
        self.Nacelle_P_ang_vel_num = create_same_shape_list(self.Nacelle_P_ang_vel)
        self.HSS_P_vel_num = create_same_shape_list(self.HSS_P_vel)
        self.HSS_P_ang_vel_num = create_same_shape_list(self.HSS_P_ang_vel)
        self.TowerTop_P_vel_num = create_same_shape_list(self.TowerTop_P_vel)
        self.TowerTop_P_ang_vel_num = create_same_shape_list(self.TowerTop_P_ang_vel)
        self.PltFm_P_vel_num = create_same_shape_list(self.PltFm_P_vel)
        self.PltFmRef_P_vel_num = create_same_shape_list(self.PltFmRef_P_vel)
        self.PltFm_P_ang_vel_num = create_same_shape_list(self.PltFm_P_ang_vel)

        for i in range(len(self.Hub_P_vel)):
            self.Hub_P_vel_num[i] = Numerify(self.Hub_P_vel[i])
            self.Hub_P_ang_vel_num[i] = Numerify(self.Hub_P_ang_vel[i])
            self.Nacelle_P_vel_num[i] = Numerify(self.Nacelle_P_vel[i])
            self.Nacelle_P_ang_vel_num[i] = Numerify(self.Nacelle_P_ang_vel[i])
            self.HSS_P_vel_num[i] = Numerify(self.HSS_P_vel[i])
            self.HSS_P_ang_vel_num[i] = Numerify(self.HSS_P_ang_vel[i])
            self.TowerTop_P_vel_num[i] = Numerify(self.TowerTop_P_vel[i])
            self.TowerTop_P_ang_vel_num[i] = Numerify(self.TowerTop_P_ang_vel[i])
            self.PltFm_P_vel_num[i] = Numerify(self.PltFm_P_vel[i])
            self.PltFm_P_ang_vel_num[i] = Numerify(self.PltFm_P_ang_vel[i])
            self.PltFmRef_P_vel_num[i]  = Numerify(self.PltFmRef_P_vel[i])

        self.TowerEle_P_vel_num = create_same_shape_list(self.TowerEle_P_vel)
        self.TowerEle_P_ang_vel_num = create_same_shape_list(self.TowerEle_P_ang_vel)

        # Tower elements----------------------------
        for i in range(len(self.TowerEle_P_vel_num)):
            for j in range(len(self.TowerEle_P_vel_num[0])):
                self.TowerEle_P_vel_num[i][j] = Numerify(self.TowerEle_P_vel[i][j])
                self.TowerEle_P_ang_vel_num[i][j] = Numerify(self.TowerEle_P_ang_vel[i][j])


def SetMultibodySystem(Parameters):
    # Blade = SetBladePara(TurbConfig, BladeData)
    # Platform = SetPlatformPara(TurbConfig)
    # TurbConfig['AllBaldeMass'] = Blade.AllBladeMass
    # Tower = SetTowerPara(TurbConfig, TowerData)
    # Nacelle = SetNacellePara(TurbConfig)
    # Drivetrain = SetDriveTrainPara(TurbConfig)
    # Hub = SetHubPara(TurbConfig)
    MBS = MBSsystem(Parameters)
    # MBS.SetDynamics(Parameters)

    return MBS


def create_same_shape_list(original_list):
    if isinstance(original_list, list):
        return [create_same_shape_list(item) for item in original_list]
    else:
        return None  # or some default value based on your requirement


def SimulateMBS(MBS, Parameters):
    dt = 0.005
    T = 10
    t = 0
    Steps = int(T / dt)

    ActiveGenSpeed = MBS.ActiveGenSpeed
    ActiveDOF = MBS.ActiveDOF

    q0 = [0, 0]
    qdt0 = [0, 0]
    qreplace = dict(zip(ActiveDOF, q0))
    ureplace = dict(zip(ActiveGenSpeed, qdt0))

    # ud_zerod = {q: 0 for q in ActiveGenSpeed}
    # q_zerod = {q: 0 for q in ActiveDOF}
    # q_zerod

    Fr_bar = []
    Frs_bar = []
    GlobalFrame = MBS.GlobalFrame
    BElmntMass = Parameters.Blade.BElmntMass
    g = 9.81

    for it in range(Steps):

        Fr = 0
        Frs = 0
        for k in range(len(ActiveDOF)):

            for i in range(Parameters.Blade.BldNodes):  # Parameters.Blade.BldNodes):
                forces = -float(BElmntMass[i] * g) * GlobalFrame.z
                masses = float(BElmntMass[i])
                torques = 0.0 * GlobalFrame.z
                inertias = 0.0

                vr = MBS.Blade1Ele_P_vel[i][k].xreplace(ureplace).xreplace(qreplace)
                Fr += vr.dot(forces)
                Rs = -masses * MBS.Blade1Ele_acc[i].xreplace(ureplace).xreplace(qreplace)
                Frs += vr.dot(Rs)

                wr = MBS.Blade1Ele_P_ang_vel[i][k].xreplace(ureplace).xreplace(qreplace)
                Fr += wr.dot(torques)
                # Ts = -(Bi.ang_acc_in(GlobalFrame).dot(inertias) +
                #       me.cross(Bi.ang_vel_in(GlobalFrame), Ii).dot(Bi.ang_vel_in(GlobalFrame)))

                vr = MBS.Blade2Ele_P_vel[i][k].xreplace(ureplace).xreplace(qreplace)
                Fr += vr.dot(forces)
                Rs = -masses * MBS.Blade2Ele_acc[i].xreplace(ureplace).xreplace(qreplace)
                Frs += vr.dot(Rs)

                wr = MBS.Blade2Ele_P_ang_vel[i][k].xreplace(ureplace).xreplace(qreplace)
                Fr += wr.dot(torques)
                # Ts = -(Bi.ang_acc_in(GlobalFrame).dot(inertias) +
                #       me.cross(Bi.ang_vel_in(GlobalFrame), Ii).dot(Bi.ang_vel_in(GlobalFrame)))

                vr = MBS.Blade3Ele_P_vel[i][k].xreplace(ureplace).xreplace(qreplace)
                Fr += vr.dot(forces)
                Rs = -masses * MBS.Blade3Ele_acc[i].xreplace(ureplace).xreplace(qreplace)
                Frs += vr.dot(Rs)

                wr = MBS.Blade3Ele_P_ang_vel[i][k].xreplace(ureplace).xreplace(qreplace)
                Fr += wr.dot(torques)
                # Ts = -(Bi.ang_acc_in(GlobalFrame).dot(inertias) +
                #       me.cross(Bi.ang_vel_in(GlobalFrame), Ii).dot(Bi.ang_vel_in(GlobalFrame)))

            Fr_bar.append(Fr)
            Frs_bar.append(Frs)

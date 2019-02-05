#!/usr/bin/env python3
import numpy as np
import math
from robot_kinematic_python.Structure import DH, Pos


class Kinematic():
    def __init__(self, Number_of_joint):
        self._total_links = Number_of_joint
        self._sDH = []
        for i in range(self._total_links):
            self._sDH.append(DH())
        self._T0 = np.eye(4)
        self._TF = np.eye(4)
        self._H0F = np.eye(4)

    # ******************************************************************
    #  * Set DH Parameters
    #  *
    #  * @index  : index
    #  * @a, d, alpha, theta0 : D-H parameter values
    #  * @active : active joints -> 1
    #  * @min    : minimum joint limit
    #  * @max    : maximum joint limit
    #  * @maxVel : maximum joint velocity limit ( >0 )
    #  *****************************************************************
    def setDH(self, index, a, d, alpha, theta0, active, min, max, maxVel):
        if index >= self._total_links:
            return -1
        else:
            self._sDH[index].alpha = alpha
            self._sDH[index].a = a
            self._sDH[index].d = d
            self._sDH[index].theta0 = theta0
            self._sDH[index].active = active
            self._sDH[index].min = min
            self._sDH[index].max = max
            self._sDH[index].maxVel = maxVel
            self._sDH[index].weight = 1.0
            self._sDH[index].theta = 0.0

            self._sDH[index].H[2, 1] = math.sin(alpha)
            self._sDH[index].H[2, 2] = math.cos(alpha)
            self._sDH[index].H[2, 3] = d
            self._sDH[index].H[3, 3] = 1.0

    # ******************************************************************
    #  * Set Base frame
    #  *****************************************************************

    def setT0(self, T):
        for i in range(4):
            for j in range(4):
                self._T0[i, j] = T[i, j]

    # ******************************************************************
    #  * Set End frame
    #  *****************************************************************
    def setTF(self, T):
        for i in range(4):
            for j in range(4):
                self._TF[i, j] = T[i, j]

    # ******************************************************************
    #  *  Make the kinematic chain is ready to use!
    #  *****************************************************************
    def readyForKinematics(self):
        self._dof = 0  # find degree of freedom
        for i in range(self._total_links):
            if self._sDH[i].active == 1:
                self._dof = self._dof+1

        #  generate active_index
        self._active_index = [0]*self._dof
        ai = 0
        for i in range(self._total_links):
            if self._sDH[i].active == 1:
                self._active_index[ai] = i
                ai = ai+1

        # find inactive link trnsformation matrix
        for i in range(self._total_links):
            if self._sDH[i].active == 1:
                c_theta = math.cos(self._sDH[i].theta0)
                s_theta = math.sin(self._sDH[i].theta0)
                c_alpha = math.cos(self._sDH[i].alpha)
                s_alpha = math.sin(self._sDH[i].alpha)

                self._sDH[i].H[0, 0] = c_theta
                self._sDH[i].H[0, 1] = -s_theta*c_alpha
                self._sDH[i].H[0, 2] = s_theta*s_alpha
                self._sDH[i].H[0, 3] = c_theta*self._sDH[i].a
                self._sDH[i].H[1, 0] = s_theta
                self._sDH[i].H[1, 1] = c_theta*c_alpha
                self._sDH[i].H[1, 2] = -c_theta*s_alpha
                self._sDH[i].H[1, 3] = s_theta*self._sDH[i].a

        self.setJoints(np.zeros((self._dof, 1)),False)

    # ******************************************************************
    #  *  Set the joitn space!
    #  *****************************************************************
    def setJoints(self, ang, check):
        i = 0
        if check==True:
            for ai in range(self._dof):
                # min max check
                i = self._active_index[ai]
                if ang[ai, 0] < self._sDH[i].min:
                    self._sDH[i].theta = self._sDH[i].min
                elif ang[ai, 0] > self._sDH[i].max:
                    self._sDH[i].theta = self._sDH[i].max
                else:
                    self._sDH[ai].theta = ang[i, 0]
        else:
            for ai in range(self._dof):
                self._sDH[ai].theta = ang[ai, 0]
        self._calFwd()

    # ******************************************************************
    #  *  Solve Forward kinematic
    #  *****************************************************************
    def _calFwd(self):
        for ai in range(self._dof):
            link_index = self._active_index[ai]
            c_theta = math.cos(self._sDH[link_index].theta +
                               self._sDH[link_index].theta0)
            s_theta = math.sin(self._sDH[link_index].theta +
                               self._sDH[link_index].theta0)
            c_alpha = math.cos(self._sDH[link_index].alpha)
            s_alpha = math.sin(self._sDH[link_index].alpha)

            self._sDH[link_index].H[0, 0] = c_theta
            self._sDH[link_index].H[0, 1] = -s_theta*c_alpha
            self._sDH[link_index].H[0, 2] = s_theta*s_alpha
            self._sDH[link_index].H[0, 3] = c_theta*self._sDH[link_index].a

            self._sDH[link_index].H[1, 0] = s_theta
            self._sDH[link_index].H[1, 1] = c_theta*c_alpha
            self._sDH[link_index].H[1, 2] = -c_theta*s_alpha
            self._sDH[link_index].H[1, 3] = s_theta*self._sDH[link_index].a
        self._sDH[0].H0i = np.dot(self._T0, self._sDH[0].H)
        for i in range(1, self._dof):
            self._sDH[i].H0i = np.dot(self._sDH[i-1].H0i, self._sDH[i].H)

        self._H0F = np.dot(self._sDH[self._total_links-1].H0i, self._TF)

    # ******************************************************************
    #  *  Get Joint position
    #  *****************************************************************
    def getJoints(self):
        ang = np.zeros((self._dof), 1)
        for ai in range(self._dof):
            ang[ai] = self._sDH[self._active_index[ai]].theta
        return ang

    # ******************************************************************
    #  *  Get End position
    #  *****************************************************************
    def getEndPos(self):
        return self._H0F[0:3, 3]

    # ******************************************************************
    #  *  Get End + Orientation
    #  *****************************************************************
    def getEndTMatrix(self):
        return self._H0F

    # ******************************************************************
    #  *  Get full jacobian
    #  *****************************************************************
    def getJacobian(self):
        i = 0
        j = 0
        ai = 0
        v1 = np.zeros(3)
        v2 = np.zeros(3)
        J = np.zeros((6, self._dof))
        for ai in range(self._dof):
            i = self._active_index[ai]
            if i == 0:
                for j in range(3):
                    v1[j] = self._T0[j, 2]
                    v2[j] = self._H0F[j, 3] - self._T0[j, 3]

                J[3, ai] = self._T0[0, 2]
                J[4, ai] = self._T0[1, 2]
                J[5, ai] = self._T0[2, 2]
            else:
                for j in range(3):
                    v1[j] = self._sDH[i-1].H0i[j, 2]
                    v2[j] = self._H0F[j, 3] - self._sDH[i-1].H0i[j, 3]

                J[3, ai] = self._sDH[i-1].H0i[0, 2]
                J[4, ai] = self._sDH[i-1].H0i[1, 2]
                J[5, ai] = self._sDH[i-1].H0i[2, 2]
            v3 = np.cross(v1, v2)
            J[0, ai] = v3[0]
            J[1, ai] = v3[1]
            J[2, ai] = v3[2]
        return J

    # ******************************************************************
    #  *  Get jacobian along position
    #  *****************************************************************
    def getJacobianPos(self):
        i = 0
        j = 0
        ai = 0
        v1 = np.zeros(3)
        v2 = np.zeros(3)
        J = np.zeros((3, self._dof))
        for ai in range(self._dof):
            i = self._active_index[ai]
            if i == 0:
                for j in range(3):
                    v1[j] = self._T0[j, 2]
                    v2[j] = self._H0F[j, 3] - self._T0[j, 3]
            else:
                for j in range(3):
                    v1[j] = self._sDH[i-1].H0i[j, 2]
                    v2[j] = self._H0F[j, 3] - self._sDH[i-1].H0i[j, 3]

            v3 = np.cross(v1, v2)
            J[0, ai] = v3[0]
            J[1, ai] = v3[1]
            J[2, ai] = v3[2]
        return J

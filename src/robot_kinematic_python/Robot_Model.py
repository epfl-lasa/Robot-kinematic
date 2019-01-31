#!/usr/bin/env python3
import numpy as np
import math
from robot_kinematic_python.Structure import DH, Pos


class Kinematic():
    def __init__(self, Number_of_joint):
        self._total_links = Number_of_joint
        self._sDH = [DH]*self._total_links
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
            return 0

    # ******************************************************************
    #  * Set Base frame
    #  *****************************************************************
    def setT0(self, T):
        for i in range(4):
            for j in range(4):
                self._T0[i, j] = T[i, j]

    # ******************************************************************
    #  *  Make the kinematic chain is ready to use!
    #  *****************************************************************
    def readyForKinematics(self):
        self._dof = 0  # find degree of freedom
        for i in range(self._total_links):
            if self._sDH[i].active == 1:
                self._dof = self._dof+1

        #  generate active_index
        self.active_index = [0]*self._dof
        ai = 0
        for i in range(self._total_links):
            if self._sDH[i].active == 1:
                ai = ai+1
                self.active_index[ai] = i

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

        self.setJoints(np.zeros((self._dof, 1)))

    # ******************************************************************
    #  *  Set the joitn space!
    #  *****************************************************************
    def setJoints(self, ang):
        i = 0
        for ai in range(self._dof):
            # min max check
            i = self.active_index[ai]
            if ang[ai] < self._sDH[i].min:
                self._sDH[i].theta = self._sDH[i].min
            elif ang[ai] > self._sDH[i].max:
                self._sDH[i].theta = self._sDH[i].max
            else:
                self._sDH[i].theta = ang[ai]
        self._calFwd()

    # ******************************************************************
    #  *  Solve Forward kinematic
    #  *****************************************************************
    def _calFwd(self):
        for ai in range(self._dof):
            link_index = self.active_index[ai]
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
            ang[ai] = self._sDH[self.active_index[ai]].theta
        return ang

    # ******************************************************************
    #  *  Get End position
    #  *****************************************************************
    def getEndPos(self):
        return self._H0F[0:3,3]

    # ******************************************************************
    #  *  Get End + Orientation
    #  *****************************************************************
    def getEndTMatrix(self):
        return self._H0F
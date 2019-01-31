#!/usr/bin/env python3
import numpy as np
import math


class Pos():
    def __init__(self):
        self.position = np.zeros((3, 1))
        self.orientation = np.zeros((3, 3))
        self.transformation = np.zeros((4, 4))


class DH():
    def __init__(self):
        self.active = 0
        self.alpha = 0.0
        self.a = 0.0
        self.d = 0.0
        self.theta0 = 0.0
        self.min = 0.0
        self.max = 0.0
        self.theta = 0.0
        self.weight = 0.0  # 0.0 <= weight <= 1.0
        self.maxVel = 0.0
        self.H = np.zeros((4, 4))
        self.H0i = np.zeros((4, 4))

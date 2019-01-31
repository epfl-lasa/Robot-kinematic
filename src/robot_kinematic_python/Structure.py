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


# Some useful constants.
M_E = 2.7182818284590452354  # e
M_LOG2E = 1.4426950408889634074  # log_2 e
M_LOG10E = 0.43429448190325182765  # log_10 e
M_LN2 = 0.69314718055994530942  # log_e 2
M_LN10 = 2.30258509299404568402  # log_e 10
M_PI = 3.14159265358979323846  # pi
M_PI_2 = 1.57079632679489661923  # pi/2
M_PI_4 = 0.78539816339744830962  # pi/4
M_1_PI = 0.31830988618379067154  # 1/pi
M_2_PI = 0.63661977236758134308  # 2/pi
M_2_SQRTPI = 1.12837916709551257390  # 2/sqrt(pi)
M_SQRT2 = 1.41421356237309504880  # sqrt(2)
M_SQRT1_2 = 0.70710678118654752440  # 1/sqrt(2)

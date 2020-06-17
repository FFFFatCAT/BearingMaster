"""
                       .::::.
                     .::::::::.
                    :::::::::::
                 ..:::::::::::'
              '::::::::::::'
                .::::::::::
           '::::::::::::::..
                ..::::::::::::.
              ``::::::::::::::::
               ::::``:::::::::'        .:::.
              ::::'   ':::::'       .::::::::.
            .::::'      ::::     .:::::::'::::.
           .:::'       :::::  .:::::::::' ':::::.
          .::'        :::::.:::::::::'      ':::::.
         .::'         ::::::::::::::'         ``::::.
     ...:::           ::::::::::::'              ``::.
    ````':.          ':::::::::'                  ::::..
                       '.:::::'                    ':'````..
"""

# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018-10-21 16:50
# @Author  : Peter.Wong
# @File    : BearingStrength.py

import os
import base64
import sys
import time
from math import cos, acos, sin, tan, asin

from xlrd import open_workbook
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox, QFileDialog
from matplotlib import pyplot as plt
from numpy import arange, array, lexsort, ones, append, savetxt
from scipy.optimize import fsolve
from scipy.special import ellipk, ellipe

import main

MainWindow_main = main.Ui_MainWindow


def cal_cp(Dw, alpha0, Dpw, ri, re, Ep, ve, tol=0.001):
    """
        calculation of the spring constant cp, acc. to ISO 16218 Function 11.
    Args:
        tol: tolerance for convergence
        Dw: diameter of the ball
        alpha0: initial contact angle
        Dpw: pitch diameter of the bearing
        ri: cross-sectional raceway groove radius, inner
        re: cross-sectional raceway groove radius, outer
        Ep: modulus of elasticity
        ve: poisson's ratio
    Returns:
        cp: float
    """
    gamma = Dw * cos(alpha0) / Dpw
    # 内外圈曲率和
    rho_i = 2 / Dw * (2 + gamma / (1 - gamma) - Dw / 2 / ri)
    rho_e = 2 / Dw * (2 - gamma / (1 + gamma) - Dw / 2 / re)
    # 内外圈曲率差
    Fip = (gamma / (1 - gamma) + Dw / 2 / ri) / (2 + gamma / (1 - gamma) - Dw / 2 / ri)
    Fep = (-gamma / (1 + gamma) + Dw / 2 / re) / (2 - gamma / (1 + gamma) - Dw / 2 / ri)
    for k in arange(0, 1, 0.001):  # 内圈迭代求解
        if k == 0:
            chi = float("inf")
        else:
            chi = 1 / k
        M = 1 - 1 / chi ** 2
        #  第一类和第二类椭圆积分
        Ki = ellipk(M)
        Ei = ellipe(M)
        Fp = 1 - 2 / (chi ** 2 - 1) * (Ki / Ei - 1)
        if abs((Fp - Fip) / Fp) < tol:
            chi_i = chi
            break
        else:
            pass
    for k in arange(0, 1, 0.001):  # 外圈迭代求解
        if k == 0:
            chi = float("inf")
        else:
            chi = 1 / k
        M = 1 - 1 / chi ** 2
        #  第一类和第二类椭圆积分
        Ke = ellipk(M)
        Ee = ellipe(M)
        Fp = 1 - 2 / (chi ** 2 - 1) * (Ke / Ee - 1)
        if abs((Fp - Fep) / Fp) < tol:
            chi_e = chi
            break
        else:
            pass
    return (
        1.48
        * Ep
        / (1 - ve ** 2)
        * (
            (
                Ki * (rho_i / chi_i ** 2 / Ei) ** (1 / 3)
                + Ke * (rho_e / chi_e ** 2 / Ee) ** (1 / 3)
            )
        )
        ** (-1.5)
    )


def cal_fc(index, alpha):
    """
        calculation of the fc acc. to ISO 281 Table 4.
        # only supports interpolation calculation at the same angle,
            the interpolation for angel between 45 and 60 need to be updated.
        # no error check functions.
    Args:
        index: float, Dw*cos(alpha)/Dpw
        alpha: float, initial contact angel.
    Returns:
        fc, float.
            if error occurred, return 'None'
    """
    arges = arange(0.01, 0.32, 0.01)
    a_45 = (
        42.1,
        51.7,
        58.2,
        63.3,
        67.3,
        70.7,
        73.5,
        75.9,
        78,
        79.7,
        81.1,
        82.3,
        83.3,
        84.1,
        84.7,
        85.1,
        85.4,
        85.5,
        85.5,
        85.4,
        85.2,
        84.9,
        84.5,
        84,
        83.4,
        82.8,
        82,
        81.3,
        80.4,
        79.6,
    )
    a_60 = (
        39.2,
        48.1,
        54.2,
        58.9,
        62.6,
        65.8,
        68.4,
        70.7,
        72.6,
        74.2,
        75.5,
        76.6,
        77.5,
        78.3,
        78.8,
        79.2,
        79.5,
        79.6,
        79.6,
        79.5,
    )
    a_75 = (37.3, 45.9, 51.7, 56.1, 59.7, 62.7, 65.2, 67.3, 69.2, 70.7)

    for i, k in enumerate(arges):
        if k <= index < arges[i + 1]:
            if alpha == acos(-1) / 4:
                return (index - k) * (a_45[i + 1] - a_45[i]) / (
                    arges[i + 1] - arges[i]
                ) + a_45[i]
            elif alpha == acos(-1) / 3:
                return (index - k) * (a_60[i + 1] - a_60[i]) / (
                    arges[i + 1] - arges[i]
                ) + a_60[i]
            elif alpha == 75 * acos(-1) / 180:
                return (index - k) * (a_75[i + 1] - a_75[i]) / (
                    arges[i + 1] - arges[i]
                ) + a_75[i]
            else:
                with open("ErrorLog.txt", "w") as f:
                    f.write(
                        time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time()))
                        + "\n"
                    )
                    f.write(
                        "Errors occurred during the calculation of coefficient fc, check please."
                    )
                    f.close()
                return None
        else:
            pass
    with open("ErrorLog.txt", "w") as f:
        f.write(time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time())) + "\n")
        f.write(
            "Errors occurred during the calculation of coefficient fc, check please."
        )
        f.close()
    return None


def cal_ca(fc, alpha, Z, Dw, bm=1.3):
    """
        calculation of the basic dynamic axial load rating according to ISO281 chapter 6.1.
    pares:
        fc: float
        alpha: float, initial contact angel
        Z: int, number of rolling elements
        Dw: float, nominal ball diameter
    return:
        Ca: float, basic dynamic axial load rating

    """
    if Dw <= 25.4:
        return [
            bm * fc * Z ** (2 / 3) * Dw ** 1.8
            if alpha == 90
            else bm * fc * Z ** (2 / 3) * Dw ** 1.8 * cos(alpha) ** 0.7 * tan(alpha)
        ][0]
    else:
        return [
            3.647 * bm * fc * Z ** (2 / 3) * Dw ** 1.4
            if alpha == acos(-1) / 2
            else 3.647
            * bm
            * fc
            * Z ** (2 / 3)
            * Dw ** 1.4
            * cos(alpha) ** 0.7
            * tan(alpha)
        ][0]


def cal_load_fre(arrays):
    """
        calculation of the equivalent load and frequency
    Args:
        arrays: [64,2], component and rad of the load
    Returns:
        array: Mx Fre_Mx -Mx Fre_-Mx N

    """
    pi = acos(-1)
    # 将数组从小（负）到大排序
    array_sorted = arrays.T[lexsort(arrays[::-1, :])].T
    load = array_sorted[0]
    rad = array_sorted[1]
    # 统计载荷值中为负的个数
    neg_pin = 0
    for i in load:
        if i < 0:
            neg_pin += 1
    n = [x / (2 * pi) for x in rad]  # 对应参考文献[1]中的'n'

    N_neg = sum(n[:neg_pin])  # 对应参考文献[1]中的'负N'
    N_pos = sum(n[neg_pin:])  # 对应参考文献[1]中的'正N'
    N = N_neg + N_pos

    sum_frequency_neg = N_neg / N  # 对应参考文献[1]中的'负frequency'
    sum_frequency_pos = N_pos / N  # 对应参考文献[1]中的'正frequency'

    load_3_n = [(load[x] ** 3) * n[x] for x in range(64)]  # 对应参考文献[1]中的'Mx^3*n'等

    sum_load_3_n_neg = sum(load_3_n[:neg_pin])  # 对应参考文献[1]中的'sum(Mx^3*n)'等
    sum_load_3_n_pos = sum(load_3_n[neg_pin:])

    eqv_load_neg = -(
        (-sum_load_3_n_neg / N_neg) ** (1 / 3)
    )  # sum_load_3_n_neg为负值，python中无法直接为负值进行开根号并返回浮点数（返回复数）
    eqv_load_pos = (sum_load_3_n_pos / N_pos) ** (1 / 3)

    return [
        eqv_load_pos,
        eqv_load_neg,
        sum_frequency_pos,
        sum_frequency_neg,
        N,
    ]


class Main(QtWidgets.QMainWindow, MainWindow_main):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        MainWindow_main.__init__(self)
        self.setupUi(self)

        self.setWindowIcon(
            QtGui.QIcon(r"D:\OneDrive\Code\BearingMaster\Data\myico.ico")
        )
        png = QtGui.QPixmap(r"D:\OneDrive\Code\BearingMaster\Data\PB.png")
        self.label_pic.setPixmap(png)

        self.cwd = os.getcwd()  # current file location
        self.fileName_choose = ""  # Lrd fatigue load
        self.pares = {}  # bearing parameters
        self.loads = {}  # extreme loads
        self.phi_ball = []  # position of the roller elements, in deg.
        self.pi = acos(-1)
        self.N = 1  # design life, total rotation in rad.
        # 创建用于储存输出文件的文件夹
        folder = self.cwd + "\\result_output"
        if not os.path.exists(folder):
            os.makedirs(folder)

        # signal
        self.buttom_ult.clicked.connect(self.cal_static)
        self.buttom_fat.clicked.connect(self.cal_life)
        self.buttom_plot.clicked.connect(self.plot_qj)
        self.buttom_read_fat_file.clicked.connect(self.open_lrd)
        # self.buttom_read_fat_file.clicked.connect(self.cal_equivalent_load)
        # self.buttom_read_FEM_file.clicked.connect(self.slot_btn_chooseFile_FEM)
        # self.buttom_fat_FEM.clicked.connect(self.cal_life_FEM)

        # timenow = time.strftime('%Y-%m-%d')
        # if timenow < '2019-05-01':
        #     pass
        # else:
        #     exit()

    def read_pares(self):
        """
            read the parameters of bearings.
            if errors occurred in reading the parameters, return {}
        Returns:
            pares: dict
        """
        self.pares = {}
        try:
            self.pares["dpw"] = float(self.line_dpw.text())
            self.pares["dc"] = float(self.line_dc.text())
            self.pares["fi"] = float(self.line_fi.text())
            self.pares["fe"] = float(self.line_fe.text())
            self.pares["z0"] = float(self.line_z0.text())
            self.pares["nz"] = float(self.line_nz.text())
            self.pares["z1"] = float(self.line_z1.text())
            self.pares["z2"] = float(self.line_z2.text())
            self.pares["alpha"] = (
                float(self.line_alpha.text()) * self.pi / 180
            )  # set unit of angle to rad.
            self.pares["ep"] = float(self.line_ep.text())
            self.pares["nu"] = float(self.line_nu.text())
            self.pares["dw"] = float(self.line_dw.text())
            self.pares["gap"] = float(self.line_gap.text())
            self.pares["gr"] = float(self.line_gr.text())
        except ValueError:
            with open("ErrorLog.txt", "a") as f:
                f.write(
                    time.strftime("%Y-%m-%d %H:%M", time.localtime(time.time())) + "\n"
                )
                f.write("Error reading bearing parameters.\n")
                f.close()
                QMessageBox.warning(
                    self,
                    "Error",
                    "<span style=' font-size:12pt;'>Error reading bearing parameters！",
                )
            return self.pares
        return self.pares

    def read_load_ult(self):
        """
            read the extreme loads Fr Fa M in KN and KNm
            if errors occurred in reading the loads, return {}.
        Returns:
            loads: dict
        """
        self.loads = {}
        try:
            self.loads["fxy"] = float(self.line_fz.text())
            self.loads["fz"] = float(self.line_fxy.text())
            self.loads["mxy"] = float(self.line_mxy.text())
        except ValueError:
            with open("ErrorLog.txt", "a") as f:
                f.write(
                    time.strftime("%Y-%m-%d %H:%M", time.localtime(time.time())) + "\n"
                )
                f.write("Errors in reading static loads.\n")
                f.close()
                QMessageBox.warning(
                    self,
                    "Error",
                    '<span style="font-size:12pt;">Error reading extreme loads！',
                )
            return self.loads
        return self.loads

    def cal_static(self):
        """
            calculation of the maximum contact force and contact pressure base on the Fr, Fa, M
            acc. to：
                "Wind Turbine Design Guideline DG03: Yaw and Pitch Rolling Bearing Life"
                "ISO 281-2007 chapter-6, Appendix B"
            ### only calculation for ball bearing supported. not support for roller bearing.

        :return:
            Qmax: float
            Smax: float
        """
        # read parameters and loads
        pares = self.read_pares()
        if not pares == {}:
            f = pares["fi"]  # 内圈滚道曲率系数
            d = pares["dw"]  # 钢球直径
            z = pares["z0"]  # 单排钢球数量
            nz = pares["nz"]  # 滚道排数
            dm = pares["dpw"]  # 钢球节圆直径
            a = pares["alpha"]  # 初始接触角
        else:
            return
        loads = self.read_load_ult()
        if not loads == {}:
            f_radial = loads["fxy"]  # 变桨轴承承受的径向力，单位为kN
            f_axial = loads["fz"]  # 变桨轴承承受的轴向力，单位为kN
            mxy = loads["mxy"]  # 变桨轴承承受的弯矩，单位为kNm
        else:
            return
        gamma = d * cos(a) / dm  # 中间变量
        # for ball bearing
        # # body a 参数
        # ra1 = 2 / d  # 较小曲率半径
        # ra2 = 2 / d  # 较大曲率半径
        # # body b 参数
        # rb1 = 2 / (d * (1 - gamma) / gamma)  # 较大曲率半径（内圈）
        # rb2 = -1 / (f * d)  # 较小曲率半径
        # r_sum = ra1 + ra2 + rb1 + rb2  # 曲率和
        # r_dif = (ra1 - ra2 + rb1 - rb2) / r_sum  # 曲率差
        # for cylindrical cross roller bearing
        # ra = 2/d
        # rb = 2/d/((1-gamma)/gamma)
        # r_sum = ra + rb
        # r_dif = 0

        r_sum = 4 / d - 1 / f / d + 2 / d * gamma / (1 - gamma)
        r_dif = (1 / f + 2 * gamma / (1 - gamma)) / (
            4 - 1 / f + 2 * gamma / (1 - gamma)
        )
        # 插值计算接触椭圆的长半轴与短半轴
        fr = [
            0,
            0.1075,
            0.3204,
            0.4795,
            0.5916,
            0.6716,
            0.7332,
            0.7948,
            0.83495,
            0.87366,
            0.90999,
            0.93657,
            0.95738,
            0.9729,
            0.983797,
            0.990902,
            0.995112,
            0.9973,
            0.9981847,
            0.9989156,
            0.9994785,
            0.9998527,
            1,
        ]
        aa = [
            1,
            1.076,
            1.2623,
            1.4556,
            1.644,
            1.8258,
            2.011,
            2.265,
            2.494,
            2.8,
            3.233,
            3.738,
            4.395,
            5.267,
            6.448,
            8.062,
            10.222,
            12.789,
            14.839,
            17.974,
            23.55,
            37.38,
            float("inf"),
        ]
        bb = [
            1,
            0.9318,
            0.8114,
            0.7278,
            0.6687,
            0.6245,
            0.5881,
            0.548,
            0.5186,
            0.4863,
            0.4499,
            0.4166,
            0.383,
            0.349,
            0.315,
            0.2814,
            0.2497,
            0.2232,
            0.2072,
            0.18822,
            0.16442,
            0.1305,
            0,
        ]

        loc = 1
        for i in range(len(fr)):
            if r_dif >= fr[i]:
                loc = i
        aa1 = (aa[loc + 1] - aa[loc]) * (r_dif - fr[loc]) / (
            fr[loc + 1] - fr[loc]
        ) + aa[loc]
        bb1 = (bb[loc + 1] - bb[loc]) * (r_dif - fr[loc]) / (
            fr[loc + 1] - fr[loc]
        ) + bb[loc]
        # 计算最大滚动体载荷与接触应力
        """
        Acc. to DG03 Appendix B.
        A 55%/45% thrust load sharing of the two rows is considered the best possible load distribution ratio because 
        of tolerances and variation of internal dimensions between the bearing rows.
        """
        if nz == 1:
            ldf = 1
        elif nz == 2:
            ldf = 0.55
        else:
            QMessageBox.warning(
                self, "Error", "<span style=' font-size:12pt;'>最多支持双列轴承计算，请检查输入参数Z0！"
            )
            return
        qmax = ldf * (
            2 * f_radial * 1000 / (z * cos(a))
            + f_axial * 1000 / (z * sin(a))
            + 4 * mxy * 1000000 / (dm * z * sin(a))
        )
        # 滚动体最大载荷
        cc = 0.0236 * aa1 * (qmax / r_sum) ** (1 / 3)  # 椭圆长半轴
        dd = 0.0236 * bb1 * (qmax / r_sum) ** (1 / 3)  # 椭圆短半轴
        smax = 1.5 * qmax / (self.pi * cc * dd)  # 滚动体最大接触应力
        self.result_qmax.setText("%.2f" % float(qmax / 1000))  # in KN
        self.result_smax.setText("%.2f" % float(smax))  # in MPa KN/mm2
        QMessageBox.information(self, "Done", '<span style=" font-size:12pt;">接触力计算完成！')
        return qmax, smax

    def cal_qj(self, pares, loads, show=True):
        """
            calculation of the contact force and contact angel for each roller element.
        Args:
            pares: {}, bearing parameters
            loads: {}, static loads
            show: bool, show the message box after calculation.
        Return:
            Q： array, contact force
            a: array, contact angel, in deg.

        """

        def fx(x):
            """
                牛顿迭代法待求解函数,通过迭代输入变量[a, r, theta]解方程，使径向力、轴向力、弯矩与输入载荷匹配。
            Args:
                x: list, [a, r, theta]
                    a: float, 外圈轴向位移
                    r: float, 外圈径向位置
                    theta: float, 角位移
            Returns:
                deviation: list, calculated [Fa, Fr, M0] according to 'x' minus the true [Fa. Fr, M0].
            """
            Fa0, Fr0, M0 = [], [], []  # initialize the [Fa, Fr, M0] for each roller
            faa, frr, ftheta = x[0], x[1], x[2]  # initialize the inputs

            for ii in arange(1, Z + 1, 1):
                fphi = self.pi * 2 / Z * ii
                fA1 = (
                    (A * sin(alpha0) + faa + Ri * ftheta * cos(fphi)) ** 2
                    + (
                        A * cos(alpha0)
                        + frr * cos(fphi)
                        - 0.5 * dc * ftheta * cos(fphi)
                    )
                    ** 2
                ) ** 0.5
                fA2 = (
                    (A * sin(alpha0) - faa - Ri * ftheta * cos(fphi)) ** 2
                    + (
                        A * cos(alpha0)
                        + frr * cos(fphi)
                        - 0.5 * dc * ftheta * cos(fphi)
                    )
                    ** 2
                ) ** 0.5
                fA3 = (
                    (A * sin(alpha0) + faa + Ri * ftheta * cos(fphi)) ** 2
                    + (
                        A * cos(alpha0)
                        + frr * cos(fphi)
                        + 0.5 * dc * ftheta * cos(fphi)
                    )
                    ** 2
                ) ** 0.5
                fA4 = (
                    (A * sin(alpha0) - faa - Ri * ftheta * cos(fphi)) ** 2
                    + (
                        A * cos(alpha0)
                        + frr * cos(fphi)
                        + 0.5 * dc * ftheta * cos(fphi)
                    )
                    ** 2
                ) ** 0.5
                # 在接触对ii处，钢球与沟道总的弹性变形
                fx1 = fA1 - A0
                fx2 = fA2 - A0
                fx3 = fA3 - A0
                fx4 = fA4 - A0
                # 内外圈发生位移后在接触对i在位置角Phi处的接触角
                sina1 = (A * sin(alpha0) + faa + Ri * ftheta * cos(fphi)) / fA1
                cosa1 = (
                    A * cos(alpha0) + frr * cos(fphi) - 0.5 * dc * ftheta * cos(fphi)
                ) / fA1
                sina2 = (A * sin(alpha0) - faa - Ri * ftheta * cos(fphi)) / fA2
                cosa2 = (
                    A * cos(alpha0) + frr * cos(fphi) - 0.5 * dc * ftheta * cos(fphi)
                ) / fA2
                sina3 = (A * sin(alpha0) + faa + Ri * ftheta * cos(fphi)) / fA3
                cosa3 = (
                    A * cos(alpha0) + frr * cos(fphi) + 0.5 * dc * ftheta * cos(fphi)
                ) / fA3
                sina4 = (A * sin(alpha0) - faa - Ri * ftheta * cos(fphi)) / fA4
                cosa4 = (
                    A * cos(alpha0) + frr * cos(fphi) + 0.5 * dc * ftheta * cos(fphi)
                ) / fA4
                # Hertz接触理论，轴向载荷与接触变形关系
                fQ1 = Kn * ((abs(fx1) + fx1) / 2) ** 1.5
                fQ2 = Kn * ((abs(fx2) + fx2) / 2) ** 1.5
                fQ3 = Kn * ((abs(fx3) + fx3) / 2) ** 1.5
                fQ4 = Kn * ((abs(fx4) + fx4) / 2) ** 1.5
                #
                Fa0.append((fQ1 * sina1 - fQ2 * sina2 + fQ3 * sina3 - fQ4 * sina4))
                Fr0.append(
                    (fQ1 * cosa1 + fQ2 * cosa2 + fQ3 * cosa3 + fQ4 * cosa4) * cos(fphi)
                )
                M0.append(
                    1
                    / 2
                    * Dpw
                    * (fQ1 * sina1 - fQ2 * sina2 + fQ3 * sina3 - fQ4 * sina4)
                    * cos(fphi)
                )
            return [Fa - sum(Fa0), Fr - sum(Fr0), Mxy - sum(M0)]

        if not pares == {}:
            alpha0 = pares["alpha"]
            Dpw = pares["dpw"]
            Dw = pares["dw"]
            dc = pares["dc"]
            fi = pares["fi"]
            fe = pares["fe"]
            ri = fi * Dw
            re = fe * Dw
            Z = pares["z0"]
            ve = pares["nu"]
            Ep = pares["ep"]
            Gr = pares["gr"]
        else:
            return  # pares is {}

        if not loads == {}:
            Fa = loads["fxy"] * 1000
            Fr = loads["fz"] * 1000
            Mxy = loads["mxy"] * 1000000
        else:
            return  # loads is {}

        # 原始沟心距
        A = (fi + fe - 1) * Dw - Gr * cos(alpha0) / 2
        A0 = (fi + fe - 1) * Dw
        Ri = Dpw / 2 + (ri - 0.5) * Dw * cos(alpha0) - Gr * (cos(alpha0)) ** 2 / 4
        # calculation of the spring constant cp
        Kn = cal_cp(Dw, alpha0, Dpw, ri, re, Ep, ve)
        # calculation of the a, r, theta
        try:
            result_fsolve = fsolve(fx, [0.1, 0.1, 0.1])
        except Exception as e:
            with open("ErrorLog.txt", "w") as f:
                f.write(
                    time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time())) + "\n"
                )
                f.write(repr(e))
                f.close()
                QMessageBox.warning(
                    self,
                    "Error",
                    "<span style=’font-size:12pt;‘>Errors occurred during calculation of the initial contact condition.",
                )
                return None
        Q11, Q22, Q33, Q44 = [], [], [], []
        a11, a22, a33, a44 = [], [], [], []
        aa, rr, theta = result_fsolve[0], result_fsolve[1], result_fsolve[2]
        self.phi_ball = []
        for i in arange(1, Z + 1, 1):
            phi = self.pi * 2 / Z * i
            self.phi_ball.append(phi / self.pi * 180)
            A1 = (
                (A * sin(alpha0) + aa + Ri * theta * cos(phi)) ** 2
                + (A * cos(alpha0) + rr * cos(phi) - 0.5 * dc * theta * cos(phi)) ** 2
            ) ** 0.5
            A2 = (
                (A * sin(alpha0) - aa - Ri * theta * cos(phi)) ** 2
                + (A * cos(alpha0) + rr * cos(phi) - 0.5 * dc * theta * cos(phi)) ** 2
            ) ** 0.5
            A3 = (
                (A * sin(alpha0) + aa + Ri * theta * cos(phi)) ** 2
                + (A * cos(alpha0) + rr * cos(phi) + 0.5 * dc * theta * cos(phi)) ** 2
            ) ** 0.5
            A4 = (
                (A * sin(alpha0) - aa - Ri * theta * cos(phi)) ** 2
                + (A * cos(alpha0) + rr * cos(phi) + 0.5 * dc * theta * cos(phi)) ** 2
            ) ** 0.5
            # 在接触对i处，钢球与沟道总的弹性变形delta_phi
            x1 = A1 - A0
            x2 = A2 - A0
            x3 = A3 - A0
            x4 = A4 - A0
            # 内外圈发生位移后在接触对i在位置角Phi处的接触角alpha_1phi
            sina1 = (A * sin(alpha0) + aa + Ri * theta * cos(phi)) / A1
            # cosa1 = (A * cos(alpha0) + rr * cos(phi) - 0.5 * dc * theta * cos(phi)) / A1
            sina2 = (A * sin(alpha0) - aa - Ri * theta * cos(phi)) / A2
            # cosa2 = (A * cos(alpha0) + rr * cos(phi) - 0.5 * dc * theta * cos(phi)) / A2
            sina3 = (A * sin(alpha0) + aa + Ri * theta * cos(phi)) / A3
            # cosa3 = (A * cos(alpha0) + rr * cos(phi) + 0.5 * dc * theta * cos(phi)) / A3
            sina4 = (A * sin(alpha0) - aa - Ri * theta * cos(phi)) / A4
            # cosa4 = (A * cos(alpha0) + rr * cos(phi) + 0.5 * dc * theta * cos(phi)) / A4
            # Hertz接触理论，轴向载荷与接触变形关系
            Q1 = Kn * ((abs(x1) + x1) / 2) ** 1.5
            Q2 = Kn * ((abs(x2) + x2) / 2) ** 1.5
            Q3 = Kn * ((abs(x3) + x3) / 2) ** 1.5
            Q4 = Kn * ((abs(x4) + x4) / 2) ** 1.5

            a11.append(asin(sina1) / self.pi * 180)
            a22.append(asin(sina2) / self.pi * 180)
            a33.append(asin(sina3) / self.pi * 180)
            a44.append(asin(sina4) / self.pi * 180)

            Q11.append(Q1)
            Q22.append(Q2)
            Q33.append(Q3)
            Q44.append(Q4)
        if show:
            QMessageBox.information(
                self, "Done", '<span style=" font-size:12pt;">接触力、接触角计算完成！'
            )
        return Q11, Q22, Q33, Q44, a11, a22, a33, a44

    def plot_qj(self):
        """
            plot the contact force vs position of the roller elements.
        Returns:
            figure
        """
        try:
            Q11, Q22, Q33, Q44, a11, a22, a33, a44 = self.cal_qj(
                self.read_pares(), self.read_load_ult()
            )
        except TypeError:  # errors in calculation of the contact force.
            return
        plt.figure()  # 创建画布
        plt.plot(
            self.phi_ball, [x / 1000 for x in Q11], color="r", label=r"Contact 1"
        )  # 将接触力转换为kN
        plt.plot(self.phi_ball, [x / 1000 for x in Q22], color="b", label=r"Contact 2")
        plt.plot(self.phi_ball, [x / 1000 for x in Q33], color="g", label=r"Contact 3")
        plt.plot(self.phi_ball, [x / 1000 for x in Q44], color="y", label=r"Contact 4")
        plt.xlabel(r"Angle / deg")
        plt.ylabel(r"Contact load / kN")
        plt.legend(loc="upper right")
        plt.show()

        plt.figure()  # 创建画布
        plt.plot(
            self.phi_ball, [x for x in a11], color="r", label=r"Contact 1"
        )  # 将接触力转换为kN
        plt.plot(self.phi_ball, [x for x in a22], color="b", label=r"Contact 2")
        plt.plot(self.phi_ball, [x for x in a33], color="g", label=r"Contact 3")
        plt.plot(self.phi_ball, [x for x in a44], color="y", label=r"Contact 4")
        plt.xlabel(r"Angle / deg")
        plt.ylabel(r"Contact angle / deg")
        plt.legend(loc="upper right")
        plt.show()

        # with open(self.cwd + "\\result_output\\Contact_load_and_Angle.txt", "w") as f:
        #     f.writelines(
        #         "Angle(deg)      Q11(N)     Q22(N)     Q33(N)     Q44(N)" + "\n"
        #     )
        #     for i in range(len(Q11)):
        #         f.writelines(
        #             "%7.1f" % (self.phi_ball[i])
        #             + "%12.1f" % (Q11[i])
        #             + "%12.1f" % (Q22[i])
        #             + "%12.1f" % (Q33[i])
        #             + "%12.1f" % (Q44[i])
        #             + "\n"
        #         )
        #     f.close()

        return

    def open_lrd(self):
        """
            open the load file.
        Returns:
            fileName_choose: string, file name and location

        """
        self.fileName_choose, file_type = QFileDialog.getOpenFileName(
            self, "选取文件", self.cwd, "excel Files (*.xlsx);;All Files (*)"  # 起始路径
        )  # 设置文件扩展名过滤,用双分号间隔
        if self.fileName_choose == "":
            return
        self.line_fat_file_dir.setText(self.fileName_choose)
        return

    # def slot_btn_chooseFile_FEM(self):
    #     # QFileDialog.getOpenFileNames得到的fileName_choose是一个列表，相对地getOpenFileName得到的直接是一个字符串
    #     fileName_choose, filetype = QFileDialog.getOpenFileNames(
    #         self, "选取文件", self.cwd, "excel Files (*.txt);;All Files (*)"  # 起始路径
    #     )  # 设置文件扩展名过滤,用双分号间隔
    #
    #     if fileName_choose == "":
    #         print("\n取消选择")
    #         return
    #
    #     print("\n你选择的文件为:")
    #     print(fileName_choose)
    #
    #     self.line_FEM_file_dir.setText(str(fileName_choose))
    #     QMessageBox.information(self, "选择的文件为", str(fileName_choose))
    #
    #     self.FEM_file_dir = fileName_choose
    #
    # def read_FEM_file(self, dir_file):
    #     """
    #
    #     :param dir_file: 传入一个文件的绝对路径，然后读取四个接触点的接触力和接触角的数据
    #     :return:
    #     """
    #     self.Q11 = []
    #     self.Q22 = []
    #     self.Q33 = []
    #     self.Q44 = []
    #     self.alpha11 = []
    #     self.alpha22 = []
    #     self.alpha33 = []
    #     self.alpha44 = []
    #     pi = acos(-1)
    #     try:
    #         file = open(dir_file, mode="r")
    #         for line in file:
    #             line_list = line.split()
    #             self.Q11.append(float(line_list[0]))
    #             self.Q22.append(float(line_list[2]))
    #             self.Q33.append(float(line_list[4]))
    #             self.Q44.append(float(line_list[6]))
    #             self.alpha11.append(float(line_list[1]) * pi / 180)
    #             self.alpha22.append(float(line_list[3]) * pi / 180)
    #             self.alpha33.append(float(line_list[5]) * pi / 180)
    #             self.alpha44.append(float(line_list[7]) * pi / 180)
    #     except Exception as e:
    #         with open("ErrorLog.txt", "w") as f:
    #             f.write(
    #                 time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time())) + "\n"
    #             )
    #             f.write(repr(e))
    #             f.close()
    #             QMessageBox.warning(
    #                 self,
    #                 "Error",
    #                 "结果文件有问题，请对应检查以下项目：\n" "数据格式共8列-- 4个接触点，每个接触点包括接触力和接触角数据，共8列",
    #             )

    def cal_equivalent_load(self):
        """
            calculation of the equivalent lrd load.
                1, read the initial lrd load cases
                2, divide the load into positive and negative parts
                3, combine all the load components
        Returns:
            combined load components and frequency, [mxy, fxy, fz, frequency]

        """
        if self.fileName_choose == "":
            self.open_lrd()
        book_Excel = open_workbook(self.fileName_choose)
        equivalent_load = {}
        try:  # read the LRD loads into array, mixed with component and frequency
            sheet = book_Excel.sheet_by_name("叶根LRD，LDD载荷")
            # excel数据读取
            list_Mx = list(map(float, sheet.col_values(0, 1, 65)))
            list_My = list(map(float, sheet.col_values(5, 1, 65)))
            list_Fx = list(map(float, sheet.col_values(15, 1, 65)))
            list_Fy = list(map(float, sheet.col_values(20, 1, 65)))
            list_Fz = list(map(float, sheet.col_values(25, 1, 65)))

            list_rad_Mx = list(map(float, sheet.col_values(0 + 2, 1, 65)))
            list_rad_My = list(map(float, sheet.col_values(5 + 2, 1, 65)))
            list_rad_Fx = list(map(float, sheet.col_values(15 + 2, 1, 65)))
            list_rad_Fy = list(map(float, sheet.col_values(20 + 2, 1, 65)))
            list_rad_Fz = list(map(float, sheet.col_values(25 + 2, 1, 65)))

            # 以数组方式输出最终提炼的数据
            equivalent_load["mx"] = array([list_Mx, list_rad_Mx])
            equivalent_load["my"] = array([list_My, list_rad_My])
            equivalent_load["fx"] = array([list_Fx, list_rad_Fx])
            equivalent_load["fy"] = array([list_Fy, list_rad_Fy])
            equivalent_load["fz"] = array([list_Fz, list_rad_Fz])
        except:
            with open("ErrorLog.txt", "w") as f:
                f.write(
                    time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time())) + "\n"
                )
                f.write("载荷数据读取出错\n")
                f.close()
                QMessageBox.warning(
                    self,
                    "Error",
                    "载荷数据页有问题，请对应检查以下项目：\n"
                    "1 载荷数据页名称为“叶根LRD，LDD载荷”;\n"
                    "2 第一行应为“Blade 1 Mx [kNm]  Time at level   Revs[rad]    Revs at level per bin[deg/s]”",
                )
            return
        # 啦啦，文件读入没有错误的时候就可以继续载荷计算啦
        eqv_load_table = ones([2, 10])  # 对应参考文献[1]中的'等效疲劳载荷'的上表
        for i, name in enumerate(["mx", "my", "fx", "fy", "fz"]):
            eqv_load = equivalent_load[name]
            eqv_load_data = cal_load_fre(eqv_load)
            eqv_load_table[0][2 * i + 0] = eqv_load_data[0]  # eqv_load_pos
            eqv_load_table[1][2 * i + 0] = eqv_load_data[1]  # eqv_load_neg
            eqv_load_table[0][2 * i + 1] = eqv_load_data[2]  # sum_frequency_pos
            eqv_load_table[1][2 * i + 1] = eqv_load_data[3]  # sum_frequency_neg
            self.N = eqv_load_data[4]  # total equivalent rotations
        # 排列组合的工况表
        load_combine = []  # 对应参考文献[1]中的'等效疲劳载荷'的下表
        for mx in range(2):
            for my in range(2):
                for fx in range(2):
                    for fy in range(2):
                        for fz in range(2):
                            load_combine = append(
                                load_combine,
                                [
                                    eqv_load_table[mx, 0],
                                    eqv_load_table[my, 2],
                                    eqv_load_table[fx, 4],
                                    eqv_load_table[fy, 6],
                                    eqv_load_table[fz, 8],
                                    eqv_load_table[mx, 1]
                                    * eqv_load_table[my, 3]
                                    * eqv_load_table[fx, 5]
                                    * eqv_load_table[fy, 7]
                                    * eqv_load_table[fz, 9],
                                ],
                            )
        load_combine = load_combine.reshape([32, 6])  # [mx, my, fx, fy, fz ,frequency]
        load_lc = ones([32, 4])
        for i in range(32):
            load_lc[i, 0] = (load_combine[i, 0] ** 2 + load_combine[i, 1] ** 2) ** (
                1 / 2
            )
            load_lc[i, 1] = (load_combine[i, 2] ** 2 + load_combine[i, 3] ** 2) ** (
                1 / 2
            )
            load_lc[i, 2] = load_combine[i, 4]
            load_lc[i, 3] = load_combine[i, 5]
        return load_lc

    def cal_l10(self, pares, loads):
        if not pares == {}:
            alpha = pares["alpha"]
            dw = pares["dw"]
            dpw = pares["dpw"]
            z = pares["z0"]
            z1 = pares["z1"]
            z2 = pares["z2"]
            ri = pares["fi"] * dw
            re = pares["fe"] * dw
            gamma = dw * cos(alpha) / dpw
        else:
            return

        try:
            fc = cal_fc(gamma, alpha)
        except TypeError:
            return
        # calculation of basic dynamic axial load rating, according to ISO281-chapter 6.1.2, only support bearing
        # with 2 raceway  now.
        Ca = (z1 + z2) * (
            (z1 / cal_ca(fc, alpha, z1, dw)) ** (10 / 3)
            + (z2 / cal_ca(fc, alpha, z2, dw)) ** (10 / 3)
        ) ** (-3 / 10)

        # calculation of basic dynamic load rating for thrust ball bearing, according to ISO16281-chapter 4.3.1.3.
        if not alpha == self.pi / 2:
            Q_ci = (
                Ca
                / (z * sin(alpha))
                * (
                    1
                    + (
                        ((1 - gamma) / (1 + gamma)) ** 1.72
                        * (ri / re * (2 * re - dw) / 2 * ri - dw) ** 0.41
                    )
                    ** (10 / 3)
                )
                ** (3 / 10)
            )
            Q_ce = (
                Ca
                / (z * sin(alpha))
                * (
                    1
                    + (
                        ((1 - gamma) / (1 + gamma)) ** 1.72
                        * (ri / re * (2 * re - dw) / 2 * ri - dw) ** 0.41
                    )
                    ** (-10 / 3)
                )
                ** (3 / 10)
            )
        else:
            Q_ci = (
                Ca
                / z
                * (1 + (((ri / re * (2 * re - dw) / 2 * ri - dw) ** 0.41) ** (10 / 3)))
                ** (3 / 10)
            )
            Q_ce = (
                Ca
                / z
                * (1 + (((ri / re * (2 * re - dw) / 2 * ri - dw) ** 0.41) ** (-10 / 3)))
                ** (3 / 10)
            )

        # calculation of dynamic equivalent rolling element load, according to ISO16281 chapter 4.3.2
        Q11, Q22, Q33, Q44, a11, a22, a33, a44 = self.cal_qj(pares, loads, show=False)
        pitch_type = self.combo_pitch_type.currentText()
        if pitch_type == "内圈变桨":
            m1 = 3
            m2 = 1 / 3
            m3 = 10 / 3
            m4 = 3 / 10
        elif pitch_type == "外圈变桨":
            m1 = 10 / 3
            m2 = 3 / 10
            m3 = 3
            m4 = 1 / 3
        Qei_R1 = (sum([Q11[i] ** m1 for i in range(len(Q11))]) / len(Q11)) ** m2
        Qeo_R1 = (sum([Q11[i] ** m3 for i in range(len(Q11))]) / len(Q11)) ** (m4)
        Qei_R2 = (sum([Q22[i] ** m1 for i in range(len(Q11))]) / len(Q11)) ** m2
        Qeo_R2 = (sum([Q22[i] ** m3 for i in range(len(Q11))]) / len(Q11)) ** (m4)
        Qei_R3 = (sum([Q33[i] ** m1 for i in range(len(Q11))]) / len(Q11)) ** m2
        Qeo_R3 = (sum([Q33[i] ** m3 for i in range(len(Q11))]) / len(Q11)) ** (m4)
        Qei_R4 = (sum([Q44[i] ** m1 for i in range(len(Q11))]) / len(Q11)) ** m2
        Qeo_R4 = (sum([Q44[i] ** 3 for i in range(len(Q11))]) / len(Q11)) ** (m4)
        # calculation of the basic dynamic load rating Qc
        if not alpha == self.pi / 2:
            Qci = (
                Ca
                / z
                / sin(alpha)
                * (
                    1
                    + (
                        ((10 - gamma) / (1 + gamma)) ** 1.72
                        * (ri / re * (2 * re - dw) / (2 * ri - dw)) ** 0.41
                    )
                    ** (10 / 3)
                )
                ** (3 / 10)
            )
            Qce = (
                Ca
                / z
                / sin(alpha)
                * (
                    1
                    + (
                        ((10 - gamma) / (1 + gamma)) ** 1.72
                        * (ri / re * (2 * re - dw) / (2 * ri - dw)) ** 0.41
                    )
                    ** (-10 / 3)
                )
                ** (3 / 10)
            )
        else:
            Qci = (
                Ca
                / z
                * (1 + ((ri / re * (2 * re - dw) / (2 * ri - dw)) ** 0.41) ** (10 / 3))
                ** (3 / 10)
            )
            Qci = (
                Ca
                / z
                * (1 + ((ri / re * (2 * re - dw) / (2 * ri - dw)) ** 0.41) ** (-10 / 3))
                ** (3 / 10)
            )
        # calculation of the basic reference rating life L10r
        L10r_R1 = ((Qci / Qei_R1) ** (-10 / 3) + (Qce / Qeo_R1) ** (-10 / 3)) ** (
            -9 / 10
        )
        L10r_R2 = ((Qci / Qei_R2) ** (-10 / 3) + (Qce / Qeo_R2) ** (-10 / 3)) ** (
            -9 / 10
        )
        L10r_R3 = ((Qci / Qei_R3) ** (-10 / 3) + (Qce / Qeo_R3) ** (-10 / 3)) ** (
            -9 / 10
        )
        L10r_R4 = ((Qci / Qei_R4) ** (-10 / 3) + (Qce / Qeo_R4) ** (-10 / 3)) ** (
            -9 / 10
        )
        # calculation of the dynamic equivalent reference load
        Pref_R1 = Ca / L10r_R1 ** (1 / 3)
        Pref_R2 = Ca / L10r_R2 ** (1 / 3)
        Pref_R3 = Ca / L10r_R3 ** (1 / 3)
        Pref_R4 = Ca / L10r_R4 ** (1 / 3)
        #
        # L10mr中间变量 , EQ(33)
        L10mr_inter_R1 = (Ca / Pref_R1) ** 3
        L10mr_inter_R2 = (Ca / Pref_R2) ** 3
        L10mr_inter_R3 = (Ca / Pref_R3) ** 3
        L10mr_inter_R4 = (Ca / Pref_R4) ** 3
        #
        # # 从计算的四个点中取出寿命最小的
        mf_1 = 0.685  # modification factor  quote:唐静--0.685是滚道硬度修正系数
        # mf_2 = 0.45  # modification factor  quote:唐静--文献看到，润滑修正系数建议0.3-0.6之间，取了中间值，和天马刚好对上，于是加上此系数0.45
        mf_2 = 1.0  # 段师傅说我们不进行修正，交给轴承厂家进行修正
        #
        L10mr = (
            min(L10mr_inter_R1, L10mr_inter_R2, L10mr_inter_R3, L10mr_inter_R4)
            * mf_1
            * mf_2
        )
        return L10mr, mf_1

    def cal_life(self):

        pares = self.read_pares()
        try:
            load_lc = self.cal_equivalent_load()
        except TypeError:
            return  # without loads for calculation, return
        list_l10mr = []
        list_l10r = []
        frequency = []

        list_l10mri = []
        list_l10ri = []

        for i in range(32):
            loads = {}
            if load_lc[i, 3] >= 0.01:
                loads_frequency = load_lc[i, 3]
                loads["mxy"] = load_lc[i, 0]
                loads["fxy"] = load_lc[i, 1]
                loads["fz"] = load_lc[i, 2]
                # Q11, Q22, Q33, Q44, a11, a22, a33, a44 = self.cal_qj(pares, loads)

                l10mri, mf_1 = self.cal_l10(pares, loads)
                l10mri_inter = loads_frequency / l10mri

                l10ri = l10mri / mf_1
                l10ri_inter = loads_frequency / l10ri

                list_l10mr.append(l10mri_inter)
                list_l10r.append(l10ri_inter)

                frequency.append(loads_frequency)
                # 每个工况的修正基本额定寿命l10mr及基本额定寿命l10r
                list_l10mri.append(l10mri)
                list_l10ri.append(l10ri)

        l10mr_weighted = 1000000 / (sum(list_l10mr))
        l10r_weighted = 1000000 / (sum(list_l10r))
        l10r_SF = l10r_weighted / self.N  # self.N 为许用寿命；safety_factor为判断准则，是否失效，小于1为失效
        l10mr_SF = (
            l10mr_weighted / self.N
        )  # self.N 为许用寿命；safety_factor为判断准则，是否失效，小于1为失效

        self.result_life_N.setText("%.0f" % float(self.N))
        self.result_l10r_life.setText("%.0f" % float(l10r_weighted))
        self.result_l10r_SF.setText("%.3f" % float(l10r_SF))
        self.result_l10mr_life.setText("%.0f" % float(l10mr_weighted))
        self.result_l10mr_SF.setText("%.3f" % float(l10mr_SF))

        # 每个工况的修正基本额定寿命l10mr及基本额定寿命l10r
        # list_l10mri_temp = [1000000 * i for i in list_l10mri]
        # list_l10ri_temp = [1000000 * i for i in list_l10ri]
        # with open(
        #         self.cwd + "\\result_output\\l10mr and l10r of each loadcases.txt", "w"
        # ) as f:
        #     f.writelines("frequency    l10r    l10mr" + "\n")
        #     for i in range(32):
        #         if self.eqv_loadcase[i, 0] >= 0.01:
        #             f.writelines(
        #                 str(round(frequency[i], 3))
        #                 + "    "
        #                 + str(int(list_l10ri_temp[i]))
        #                 + "    "
        #                 + str(int(list_l10mri_temp[i]))
        #                 + "\n"
        #             )
        #
        #     f.writelines("\n" + "------------总体寿命--------------" + "\n")
        #     f.writelines("许用寿命 N     " + "%.0f" % float(self.N) + "\n")
        #     f.writelines("L10r寿命       " + "%.0f" % float(l10r_weighted) + "\n")
        #     f.writelines("L10r安全系数   " + "%.3f" % float(l10r_SF) + "\n")
        #     f.writelines("L10mr寿命      " + "%.0f" % float(l10mr_weighted) + "\n")
        #     f.writelines("L10mr安全系数   " + "%.3f" % float(l10mr_SF) + "\n")
        #     f.close()
        return

    def cal_l10_FEM(self):
        pi = acos(-1)
        alpha0 = self.alpha * pi / 180
        Dpw = self.dpw
        Dw = self.dw
        dc = self.dc
        gamma = Dw * cos(alpha0) / Dpw
        bm = 1.3
        fc = 51.7
        fi = self.fi
        fe = self.fe
        ri = fi * Dw
        re = fe * Dw
        Z = self.z0
        Z1 = self.z1
        Z2 = self.z2
        ii = 2
        ve = self.nu
        Ep = self.ep

        # 计算当量动载荷 - 参考ISO16281
        Qei_R1 = (
            sum([self.Q11[i] ** (10 / 3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (3 / 10)
        Qeo_R1 = (
            sum([self.Q11[i] ** (3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (1 / 3)
        Qei_R2 = (
            sum([self.Q22[i] ** (10 / 3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (3 / 10)
        Qeo_R2 = (
            sum([self.Q22[i] ** (3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (1 / 3)
        Qei_R3 = (
            sum([self.Q33[i] ** (10 / 3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (3 / 10)
        Qeo_R3 = (
            sum([self.Q33[i] ** (3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (1 / 3)
        Qei_R4 = (
            sum([self.Q44[i] ** (10 / 3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (3 / 10)
        Qeo_R4 = (
            sum([self.Q44[i] ** (3) for i in range(len(self.Q11))]) / len(self.Q11)
        ) ** (1 / 3)
        # 计算基本额定动载荷 - 参考ISO16281
        Ca1 = (
            3.647
            * fc
            * bm
            * (cos(alpha0)) ** 0.7
            * Z ** (2 / 3)
            * Dw ** (1.4)
            * tan(alpha0)
        )
        Ca2 = Ca1
        Ca1_Z1 = (Z1 / Ca1) ** (10 / 3)
        Ca2_Z2 = (Z2 / Ca2) ** (10 / 3)
        Ca_Z_1_2 = (Ca1_Z1) + (Ca2_Z2)
        Ca_Z = (Ca_Z_1_2) ** (-3 / 10)
        Z11 = Z1 + Z2
        Ca = (Ca_Z) * Z11
        m1_i = ((1 - gamma) / (1 + gamma)) ** 1.72 * (
            ri / re * ((2 * re - Dw) / (2 * ri - Dw)) ** 0.41
        )
        m2_i = (1 + m1_i ** (10 / 3)) ** (3 / 10)
        m2_ii = (m2_i) / (Z * sin(alpha0))
        Qci = Ca * (m2_ii)
        m1_e = ((1 - gamma) / (1 + gamma)) ** 1.72 * (
            ri / re * ((2 * re - Dw) / (2 * ri - Dw)) ** 0.41
        )
        m2_e = (1 + m1_e ** (-10 / 3)) ** (3 / 10)
        m2_ee = (m2_e) / (Z * sin(alpha0))
        Qce = Ca * (m2_ee)
        # L10r计算 - 参考ISO16281,EQ(29)
        L10r_R1 = ((Qci / Qei_R1) ** (-10 / 3) + (Qce / Qeo_R1) ** (-10 / 3)) ** (
            -9 / 10
        )
        L10r_R2 = ((Qci / Qei_R2) ** (-10 / 3) + (Qce / Qeo_R2) ** (-10 / 3)) ** (
            -9 / 10
        )
        L10r_R3 = ((Qci / Qei_R3) ** (-10 / 3) + (Qce / Qeo_R3) ** (-10 / 3)) ** (
            -9 / 10
        )
        L10r_R4 = ((Qci / Qei_R4) ** (-10 / 3) + (Qce / Qeo_R4) ** (-10 / 3)) ** (
            -9 / 10
        )
        # L10mr计算 - 参考ISO16281,EQ(31)
        Pref_R1 = Ca / (L10r_R1) ** (1 / 3)
        Pref_R2 = Ca / (L10r_R2) ** (1 / 3)
        Pref_R3 = Ca / (L10r_R3) ** (1 / 3)
        Pref_R4 = Ca / (L10r_R4) ** (1 / 3)

        # L10mr中间变量 , EQ(33)
        L10mr_inter_R1 = (Ca / Pref_R1) ** 3
        L10mr_inter_R2 = (Ca / Pref_R2) ** 3
        L10mr_inter_R3 = (Ca / Pref_R3) ** 3
        L10mr_inter_R4 = (Ca / Pref_R4) ** 3

        # 从计算的四个点中取出寿命最小的
        mf_1 = 0.685  # modification factor  quote:唐静--0.685是滚道硬度修正系数好像
        mf_2 = 0.45  # modification factor  quote:唐静--文献看到，润滑修正系数建议0.3-0.6之间，取了中间值，和天马刚好对上，于是加上此系数

        L10mr = (
            min(L10mr_inter_R1, L10mr_inter_R2, L10mr_inter_R3, L10mr_inter_R4)
            * mf_1
            * mf_2
        )
        # self.result_l10.setText("%.3f" % float(L10mr))
        # if ult_or_fat == 0:
        #     QMessageBox.information(self, 'Done',
        #                             '<span style=\" font-size:12pt;\">接触力计算完成！')
        # else:
        #     pass
        #
        # QMessageBox.information(self, 'Done',
        #                         '<span style=\" font-size:12pt;\">疲劳寿命计算完成！')
        return [L10mr, mf_1]

    def cal_life_FEM(self):
        list_l10mr = []
        list_l10r = []
        frequency = []

        list_l10mri = []
        list_l10ri = []

        FEM_num = len(self.FEM_file_dir)

        _ = self.fat_read()

        try:
            for i in range(FEM_num):
                loadcase_frequency = self.eqv_loadcase[i, 0]
                self.mxy = self.eqv_loadcase[i, 3]
                self.fxy = self.eqv_loadcase[i, 6]
                self.fz = self.eqv_loadcase[i, 7]
                # _ = self.cal_qj_FEM()
                _ = self.read_FEM_file(self.FEM_file_dir[i])

                l10mri = self.cal_l10_FEM()[0]
                mf_1 = self.cal_l10_FEM()[1]
                l10mri_inter = loadcase_frequency / l10mri

                l10ri = l10mri / mf_1
                l10ri_inter = loadcase_frequency / l10ri

                list_l10mr.append(l10mri_inter)
                list_l10r.append(l10ri_inter)
                frequency.append(loadcase_frequency)
                # 每个工况的修正基本额定寿命l10mr及基本额定寿命l10r
                list_l10mri.append(l10mri)
                list_l10ri.append(l10ri)

            l10mr_weighted = 1000000 * sum(frequency) / (sum(list_l10mr))
            l10r_weighted = 1000000 * sum(frequency) / (sum(list_l10r))
            l10r_SF = (
                l10r_weighted / self.N
            )  # self.N 为许用寿命；safety_factor为判断准则，是否失效，小于1为失效
            l10mr_SF = (
                l10mr_weighted / self.N
            )  # self.N 为许用寿命；safety_factor为判断准则，是否失效，小于1为失效

            self.result_life_N.setText("%.0f" % float(self.N))
            self.result_l10r_life.setText("%.0f" % float(l10r_weighted))
            self.result_l10r_SF.setText("%.3f" % float(l10r_SF))
            self.result_l10mr_life.setText("%.0f" % float(l10mr_weighted))
            self.result_l10mr_SF.setText("%.3f" % float(l10mr_SF))

            # 每个工况的修正基本额定寿命l10mr及基本额定寿命l10r
            list_l10mri_temp = [1000000 * i for i in list_l10mri]
            list_l10ri_temp = [1000000 * i for i in list_l10ri]
            with open(
                self.cwd + "\\result_output\\l10mr and l10r of each loadcases_FEM.txt",
                "w",
            ) as f:
                f.writelines("frequency    l10r    l10mr    result_file_name" + "\n")
                for i in range(FEM_num):
                    f.writelines(
                        str(round(frequency[i], 3))
                        + "    "
                        + str(int(list_l10ri_temp[i]))
                        + "    "
                        + str(int(list_l10mri_temp[i]))
                        + "    "
                        + str(self.FEM_file_dir[i])
                        + "\n"
                    )

                f.writelines("\n" + "------------总体寿命--------------" + "\n")
                f.writelines("许用寿命 N     " + "%.0f" % float(self.N) + "\n")
                f.writelines("L10r寿命       " + "%.0f" % float(l10r_weighted) + "\n")
                f.writelines("L10r安全系数   " + "%.3f" % float(l10r_SF) + "\n")
                f.writelines("L10mr寿命      " + "%.0f" % float(l10mr_weighted) + "\n")
                f.writelines("L10mr安全系数   " + "%.3f" % float(l10mr_SF) + "\n")
                f.close()

        except Exception as e:
            with open("ErrorLog.txt", "w") as f:
                f.write(
                    time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time())) + "\n"
                )
                f.write(repr(e))
                f.close()
                QMessageBox.warning(self, "Error", "请先读取载荷文件！")

        return


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Main()  # 创建QT对象
    window.show()  # QT对象显示
    sys.exit(app.exec_())

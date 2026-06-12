from manim import *
from scipy.spatial.transform import Rotation
import numpy as np
import simutils as su

data = su.extract()
t = data["data"].t
x, y, z = data["data"].x, data["data"].y, data["data"].z
Ox, Oy, Oz = data["data"].Ox, data["data"].Oy, data["data"].Oz
q0, q1, q2, q3 = data["data"].q0, data["data"].q1, data["data"].q2, data["data"].q3
T_trans, T_rot, U_gr, U_em, E = data["data"].T_trans, data["data"].T_rot, data["data"].U_gr, data["data"].U_em, data["data"].E
mu_x, mu_y, mu_z = data["data"].mu_x, data["data"].mu_y, data["data"].mu_z
SCALE = 100
SLOPE = np.tan(np.deg2rad(10))
RADIUS = 0.005
NORMAL = np.array([0,-SLOPE, 1]) / np.linalg.norm(np.array([0, -SLOPE, 1]))

def quat_mat(Q0,Q1,Q2,Q3):
    return Rotation.from_quat([Q1,Q2,Q3,Q0]).as_matrix()


def mu_rot(t):
    mu_body = np.array([mu_x[t], mu_y[t], mu_z[t]])
    rot = quat_mat(q0[t], q1[t], q2[t], q3[t])
    return rot @ mu_body


def calc_com(t):
    pos = np.array([x[t], y[t], SLOPE * y[t]])
    return pos + RADIUS * NORMAL


class Step1(ThreeDScene):
    def construct(self):
        surface = Surface(
            lambda u, v: SCALE * np.array([u, v, v * np.tan(np.deg2rad(10))]),
            u_range=[-0.04, 0.04],
            v_range=[-0.04, 0.04],
            resolution=8,
            fill_opacity=0.5,
        )
        time = ValueTracker(0)
        ball = Sphere(
            radius = SCALE * RADIUS,
            resolution = (50, 50)
        ).set_color(RED)
        def update_ball(mobject):
            t = time.get_value()
            pos = calc_com(t)
            mobject.move_to(SCALE * pos)
        ball.add_updater(update_ball)
        mu = Arrow3D(
            start = ball.get_center(),
            end = ball.get_center() + SCALE * mu_rot(0),
            resolution = 20
        )
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=-45 * DEGREES,
            frame_center=ORIGIN,
        )
        #self.camera.set_zoom(2.0)
        self.add(surface, ball, mu)
        self.play(time.animate.set_value(len(x)-1), run_time=10, rate_func=linear)
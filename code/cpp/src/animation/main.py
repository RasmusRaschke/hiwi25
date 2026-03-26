from manim import *
import numpy as np

phi = np.deg2rad(15)
from pathlib import Path
from types import SimpleNamespace

def load(path):
    """Extracts data with headers from a csv file

    Parameters
    ----------
    path : string
        Absolute or relative path to csv; bash-style wildcards are possible
    """
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    headers = data.dtype.names
    names = [h.strip() for h in headers]
    cols = {clean: data[orig] for clean, orig in zip(names, headers)}
    return cols
folder = Path(".")
datasets = {}
for csv in folder.glob("*.csv"):
    key = csv.stem
    cols = load(csv)
    datasets[key] = SimpleNamespace(**cols)

t = datasets["ramp"].t; x = datasets["ramp"].x; y = datasets["ramp"].y; z = datasets["ramp"].z
vx = datasets["ramp"].vx; vy = datasets["ramp"].vy; vz = datasets["ramp"].vz
qw = datasets["ramp"].q0; qx = datasets["ramp"].q1; qy = datasets["ramp"].q2; qz = datasets["ramp"].q3
quat = np.column_stack([qw, qx, qy, qz])
quat /= np.linalg.norm(quat, axis=1, keepdims=True)
B = np.array([0.0, 0.0, 1.0])
mu = np.array([0.0, 1.0, 0.0])
B_scale = 0.9; mu_scale = 0.9; vel_scale = 0.35

def quat_to_rot(q):
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([
        [1 - 2*(y*y + z*z),     2*(x*y - z*w),     2*(x*z + y*w)],
        [    2*(x*y + z*w), 1 - 2*(x*x + z*z),     2*(y*z - x*w)],
        [    2*(x*z - y*w),     2*(y*z + x*w), 1 - 2*(x*x + y*y)]
    ])

def lin_interpol(q0, q1, eps):
    q0 /= np.linalg.norm(q0)
    q1 /= np.linalg.norm(q1)
    scal = np.dot(q0, q1)
    if scal < 0.0:
        q1 = -q1
        scal = -scal
    scal = np.clip(scal, -1.0, 1.0)
    if scal > 0.9995:
        q = (1 - eps) * q0 + eps * q1
        return q / np.linalg.norm(q)
    
    theta0 = np.arccos(scal)
    sin_theta0 = np.sin(theta0)
    theta = theta0 * eps
    s0 = np.sin(theta0 - theta) / sin_theta0
    s1 = np.sin(theta) / sin_theta0

    return s0 * q0 + s1 * q1

def interpol_scal(tt, arr):
    return np.interp(tt, t, arr)

def interpol_quat(tt):
    if tt <= t[0]:
        return quat[0]
    if tt >= t[-1]:
        return quat[-1]

    i = np.searchsorted(t, tt) -1
    i = np.clip(i, 0, len(t)-2)
    t0, t1 = t[i], t[i+1]
    eps = (tt - t0) / (t1 - t0)
    return lin_interpol(quat[i], quat[i+1], eps)

def state(tt):
    pos = np.array([interpol_scal(tt, x), interpol_scal(tt, y), interpol_scal(tt, z)], dtype=float)
    vel = np.array([interpol_scal(tt, vx), interpol_scal(tt, vy), interpol_scal(tt, vz)], dtype=float)
    q = interpol_quat(tt)
    return pos, vel, q

class MagneticSphere(ThreeDScene):
    def construct(self):
        resolution_fa = 24
        self.set_camera_orientation(phi=75 * DEGREES, theta=330 * DEGREES)

        def param_surface(u, v):
            x = u
            y = v
            z = - np.tan(phi) * y
            return np.array([x, y, z])
        
        plane = Surface(
            param_surface,
            resolution=(resolution_fa, resolution_fa),
            v_range=[-2, 2],
            u_range=[-2, 2]
        )

        plane.scale(2, about_point=ORIGIN)
        plane.set_style(fill_opacity=1, stroke_color=GREEN)
        plane.set_fill(BLUE, opacity=0.5)

        sphere = Surface(
            lambda u, v: np.array([
                1.5 * np.cos(u) * np.cos(v),
                1.5 * np.cos(u) * np.sin(v),
                1.5 * np.sin(u)
            ]), v_range=[0, TAU], u_range=[-PI /2, PI/2],
            checkerboard_colors=[RED_D, RED_E], resolution=(15,32)
        )
        self.renderer.camera.light_source.move_to(3*IN)
        time_tracker = ValueTracker(t[0])
        sphere.save_state()
        def update_sphere(mob):
            tt = time_tracker.get_value()
            pos, vel, q = state(tt)
            R = quat_to_rot(q)
            mob.restore()
            mob.apply_matrix(R)
            mob.shift(pos)
        sphere.add_updater(update_sphere)
        com = Dot3D(radius=0.04, color=WHITE)
        com.add_updater(
            lambda m: m.move_to(state(time_tracker.get_value())[0])
        )

        trajectory = TracedPath(
            com.get_center,
            stroke_color=WHITE,
            stroke_width=4,
            dissipating_time=0,
        )

        vel_vector = always_redraw(
            lambda: (
                lambda pos, vel: Arrow3D(pos, pos + vel_scale * vel, color=YELLOW, thickness=0.02,)
            )(*state(time_tracker.get_value())[:2])
        )

        B_vector = always_redraw(
            lambda: (
                lambda pos: Arrow3D(pos, pos+ B_scale * B, color=BLUE, thickness=0.02,)
            )(state(time_tracker.get_value())[0])
        )

        mu_vector = always_redraw(
            lambda: (
                lambda pos, q: Arrow3D(pos, pos+ mu_scale * (quat_to_rot(q) @ mu), color=RED, thickness=0.02,)
            )(
                state(time_tracker.get_value())[0],
                state(time_tracker.get_value())[2],
            )
        )

        axes = ThreeDAxes()
        self.add(plane, sphere, vel_vector, com, B_vector, mu_vector)
        self.play(
            time_tracker.animate.set_value(t[-1]),
            run_time=t[-1] - t[0],
            rate_func = linear,
        )
        self.wait(10)
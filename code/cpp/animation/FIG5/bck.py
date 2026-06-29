from manim import *
from scipy.spatial.transform import Rotation
from scipy.interpolate import PchipInterpolator
import numpy as np
import math as m
import simutils as su

data = su.extract()
t = data["australia"].t

datasets = [
    {"name": "C", "x": data["australia"].x, "y": data["australia"].y, "z": data["australia"].z,
     "mu_x": data["australia"].mu_x, "mu_y": data["australia"].mu_y, "mu_z": data["australia"].mu_z, "color": RED},
    {"name": "H", "x": data["hamburg"].x, "y": data["hamburg"].y, "z": data["hamburg"].z,
     "mu_x": data["hamburg"].mu_x, "mu_y": data["hamburg"].mu_y, "mu_z": data["hamburg"].mu_z, "color": BLUE},
    {"name": "J", "x": data["jakarta"].x, "y": data["jakarta"].y, "z": data["jakarta"].z,
     "mu_x": data["jakarta"].mu_x, "mu_y": data["jakarta"].mu_y, "mu_z": data["jakarta"].mu_z, "color": YELLOW},
    {"name": "T", "x": data["tokyo"].x, "y": data["tokyo"].y, "z": data["tokyo"].z,
     "mu_x": data["tokyo"].mu_x, "mu_y": data["tokyo"].mu_y, "mu_z": data["tokyo"].mu_z, "color": GREEN},
]

print(datasets[0]["y"][-1])

SCALE = 50
SLOPE = np.tan(np.deg2rad(20))
RADIUS = 0.005
NORMAL = np.array([0, -SLOPE, 1]) / np.linalg.norm([0, -SLOPE, 1])
mu_display_scale = 0.02

for d in datasets:
    d["x_s"] = PchipInterpolator(t, d["x"])
    d["y_s"] = PchipInterpolator(t, d["y"])
    d["z_s"] = PchipInterpolator(t, d["z"])
    d["mx_s"] = PchipInterpolator(t, d["mu_x"])
    d["my_s"] = PchipInterpolator(t, d["mu_y"])
    d["mz_s"] = PchipInterpolator(t, d["mu_z"])

def quat_mat(q0, q1, q2, q3):
    return Rotation.from_quat([q1, q2, q3, q0]).as_matrix()

def state_at(sim_time, d):
    tt = np.clip(sim_time, t[0], t[-1])

    x = d["x_s"](tt)
    y = d["y_s"](tt)
    pos = np.array([x, y, SLOPE * y]) + RADIUS * NORMAL

    mu_body = np.array([d["mx_s"](tt), d["my_s"](tt), d["mz_s"](tt)])
    return pos, mu_body

class Step1(ThreeDScene):
    def construct(self):
        time = ValueTracker(t[0])

        # choose the static frame used to define the scene center
        frame_t = t[round(3 * (len(t) / 4))]
        all_points = np.array([state_at(frame_t, d)[0] for d in datasets])
        scene_center = SCALE * all_points.mean(axis=0)

        def world_pos(com):
            return SCALE * com - scene_center

        self.set_camera_orientation(
            phi=60 * DEGREES,
            theta=-45 * DEGREES,
            zoom=0.20,
        )

        balls = VGroup()
        labels = VGroup()
        trails = VGroup()
        mu_arrows = VGroup()
        mu_trails = VGroup()

        for d in datasets:
            ball = Sphere(radius=SCALE * RADIUS, resolution=(50, 50), fill_opacity=1.0).set_color(d["color"])
            ball.add_updater(lambda m, d=d: m.move_to(world_pos(state_at(time.get_value(), d)[0])))

            label = always_redraw(
                lambda d=d, b=ball: Text(d["name"], font_size=10).next_to(b, LEFT, buff=0.04)
            )
            self.add_fixed_orientation_mobjects(label)

            trail = TracedPath(
                ball.get_center,
                stroke_width=2,
                stroke_color=d["color"],
                dissipating_time=None,
            )

            def make_mu_arrow(d=d):
                com, mu_body = state_at(time.get_value(), d)
                n = np.linalg.norm(mu_body)
                if n < 1e-8:
                    mu_body = np.array([0, 0, 1.0])
                    n = 1.0
                start = world_pos(com)
                end = start + SCALE * mu_display_scale * mu_body / n
                return Arrow3D(
                    start=start,
                    end=end,
                    resolution=16,
                    color=d["color"],
                )

            mu_arrow = always_redraw(make_mu_arrow)

            mu_trail = TracedPath(
                mu_arrow.get_end,
                stroke_width=1,
                stroke_color=d["color"],
                dissipating_time=None,
            )

            balls.add(ball)
            labels.add(label)
            trails.add(trail)
            mu_arrows.add(mu_arrow)
            mu_trails.add(mu_trail)

        surface = Surface(
            lambda u, v: np.array([u, v, v * np.tan(np.deg2rad(10))]) * SCALE - scene_center,
            u_range=[-0.08, 0.08],
            v_range=[-0.6, 0.01],
            resolution=15,
            fill_opacity=0.5,
        )

        self.add(surface, balls, labels, trails, mu_arrows, mu_trails)

        # animate after centering has been fixed
        self.play(time.animate.set_value(t[-1]), run_time=10, rate_func=linear)
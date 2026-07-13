from manim import *
from scipy.spatial.transform import Rotation
from scipy.interpolate import PchipInterpolator
import numpy as np
import simutils as su

# ----------------------------
# Data
# ----------------------------
data = su.extract()
t = data["sample1"].t

datasets = [
    {
        "name": "S1",
        "x": data["sample1"].x, "y": data["sample1"].y, "z": data["sample1"].z,
        "mu_x": data["sample1"].mu_x, "mu_y": data["sample1"].mu_y, "mu_z": data["sample1"].mu_z,
        "color": BLUE,
    },
    {
        "name": "S2",
        "x": data["sample2"].x, "y": data["sample2"].y, "z": data["sample2"].z,
        "mu_x": data["sample2"].mu_x, "mu_y": data["sample2"].mu_y, "mu_z": data["sample2"].mu_z,
        "color": YELLOW,
    },
    {
        "name": "S3",
        "x": data["sample3"].x, "y": data["sample3"].y, "z": data["sample3"].z,
        "mu_x": data["sample3"].mu_x, "mu_y": data["sample3"].mu_y, "mu_z": data["sample3"].mu_z,
        "color": RED,
    },
    {
        "name": "S4",
        "x": data["sample4"].x, "y": data["sample4"].y, "z": data["sample4"].z,
        "mu_x": data["sample4"].mu_x, "mu_y": data["sample4"].mu_y, "mu_z": data["sample4"].mu_z,
        "color": PINK,
    },
    
]

SCALE = 50
SLOPE = np.tan(np.deg2rad(10))
RADIUS = 0.005
NORMAL = np.array([0, -SLOPE, 1]) / np.linalg.norm([0, -SLOPE, 1])
mu_display_scale = 0.02

# ----------------------------
# Parameters
# ----------------------------
EXPORT_STILL = False
EXPORT_TIME = 3.7 * t[-1] / 4
MU_TRAIL_DURATION = 1.5
SHOW_TIME_TRACKER = False
ZOOM = 0.25
ANIMATION_END = 0.65  # seconds
# ----------------------------
# Interpolators
# ----------------------------
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

def hud_point(scene, x, y):
    R = scene.camera.get_rotation_matrix()
    return scene.camera.frame_center + np.linalg.inv(R) @ np.array([x, y, 0.0])

class Anim(ThreeDScene):
    def construct(self):
        time = ValueTracker(EXPORT_TIME if EXPORT_STILL else t[0])

        # Choose one reference frame so the whole scene is centered consistently.
        frame_t = EXPORT_TIME if EXPORT_STILL else t[int(0.75 * (len(t) - 1))]
        all_points = np.array([state_at(frame_t, d)[0] for d in datasets])
        scene_center = SCALE * all_points.mean(axis=0)

        def world_pos(com):
            return SCALE * com - scene_center

        self.set_camera_orientation(
            phi=65 * DEGREES,
            theta=-45 * DEGREES,
            zoom=ZOOM,
        )

        balls = VGroup()
        labels = VGroup()
        trails = VGroup()
        mu_arrows = VGroup()
        mu_trails = VGroup()

        surface = Surface(
            lambda u, v: np.array([u, v, v * SLOPE]) * SCALE - scene_center,
            u_range=[-0.12, 0.12],
            v_range=[-0.6, 0.01],
            resolution=15,
            fill_opacity=0.5,
        )

        def static_trail(d, t_end, n_points=200):
            ts = np.linspace(t[0], t_end, n_points)
            pts = [world_pos(state_at(tt, d)[0]) for tt in ts]
            path = VMobject()
            if len(pts) >= 2:
                path.set_points_smoothly(pts)
            else:
                path.set_points_as_corners(pts)
            path.set_stroke(color=d["color"], width=2)
            return path

        for d in datasets:
            # Ball
            ball = Sphere(
                radius=SCALE * RADIUS,
                resolution=(50, 50),
                fill_opacity=1.0,
                stroke_width=0,
            ).set_color(d["color"])

            if EXPORT_STILL:
                com, _ = state_at(EXPORT_TIME, d)
                ball.move_to(world_pos(com))
            else:
                ball.add_updater(
                    lambda m, d=d: m.move_to(world_pos(state_at(time.get_value(), d)[0]))
                )

            # Label
            if EXPORT_STILL:
                label = Text(d["name"], font_size=10, color=BLACK)
                label.next_to(ball, LEFT, buff=0.30)
                label.add_background_rectangle(color=WHITE, opacity=0.75, buff=0.01)
            else:
                label = always_redraw(
                    lambda d=d, b=ball: (
                        Text(d["name"], font_size=5, color=BLACK)
                        .next_to(b, LEFT, buff=0.30)
                        .add_background_rectangle(color=WHITE, opacity=0.5, buff=0.01)
                    )
                )

            self.add_fixed_orientation_mobjects(label)
            self.add_foreground_mobjects(label)


            # Mu arrow
            def make_mu_arrow(d=d):
                com, mu_body = state_at(time.get_value(), d)
                n = np.linalg.norm(mu_body)
                if n < 1e-8:
                    mu_body = np.array([0, 0, 1.0])
                    n = 1.0
                start = world_pos(com)
                end = start + SCALE * mu_display_scale * mu_body / n
                return Arrow3D(start=start, end=end, resolution=30, color=d["color"])

            mu_arrow = always_redraw(make_mu_arrow)
            # Trail
            if EXPORT_STILL:
                trail = static_trail(d, EXPORT_TIME)
            else:
                trail = TracedPath(
                    ball.get_center,
                    stroke_width=2,
                    stroke_color=d["color"],
                    dissipating_time=None,
                )

                mu_trail = TracedPath(
                    mu_arrow.get_end,
                    stroke_width=1,
                    stroke_color=d["color"],
                    dissipating_time=MU_TRAIL_DURATION,
                )
                mu_trails.add(mu_trail)
            
            balls.add(ball)
            #labels.add(label)
            trails.add(trail)
            mu_arrows.add(mu_arrow)

        self.add(surface, balls, mu_arrows, trails)

        if not EXPORT_STILL:
            self.add(mu_trails)

        """if SHOW_TIME_TRACKER:
            time_value = DecimalNumber(
                EXPORT_TIME if EXPORT_STILL else time.get_value(),
                num_decimal_places=2,
                font_size=32,
                color=BLACK,
            )
            if not EXPORT_STILL:
                time_value.add_updater(lambda m: m.set_value(time.get_value()))

            time_text = Text("t =", font_size=32, color=BLACK)
            hud_text = VGroup(time_text, time_value).arrange(RIGHT, buff=0.15)

            box = SurroundingRectangle(
                hud_text,
                buff=0.12,
                color=BLACK,
                fill_color=WHITE,
                fill_opacity=0.9,
                stroke_width=1,
            )

            hud = VGroup(box, hud_text).to_corner(UL, buff=0.25)

            self.add_fixed_in_frame_mobjects(hud)"""
        if EXPORT_STILL:
            self.wait(0.1)
        else:
            self.play(time.animate.set_value(ANIMATION_END), run_time=15, rate_func=linear)
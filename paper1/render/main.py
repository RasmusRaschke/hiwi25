from manim import *
from scipy.spatial.transform import Rotation
from scipy.interpolate import PchipInterpolator
import numpy as np
import simutils as su

# ----------------------------
# Data
# Use this dictionary to import all necessary data as indicated. The color will be mapped to the color of the sphere.  
# ----------------------------
data = su.extract()
t = data["tornado"].t

datasets = [
    {
        "name": "S",
        "x": data["tornado"].x, "y": data["tornado"].y, "z": data["tornado"].z,
        "mu_x": data["tornado"].mu_x, "mu_y": data["tornado"].mu_y, "mu_z": data["tornado"].mu_z,
        "color": RED,
    },
]

# ----------------------------
# Parameters
# ----------------------------
SCALE = 50 # Do not change this, rather change zoom.
SLOPE = np.tan(np.deg2rad(5)) # Add slope (opening angle) of the incline in degree.
# slope boundaries:
SLOPE_X_MIN = -0.05
SLOPE_X_MAX = 0.05
SLOPE_Y_MIN = -0.05
SLOPE_Y_MAX = 0.05
SURFACE_OP = 0.5 # Slope opacity
RADIUS = 0.005 # Radius of the ball in cm
SPHERE_OP = 1.0 # Sphere opacity
NORMAL = np.array([0, -SLOPE, 1]) / np.linalg.norm([0, -SLOPE, 1]) # Do not change
mu_display_scale = 0.02 # Determines how long the magnetic moment vector is
EXPORT_STILL = False # If TRUE: Exports a single freeze frame at time=EXPORT_TIME
EXPORT_TIME = 3.7 * t[-1] / 4 # see above
MU_TRAIL_DURATION = 1.5 # Duration after which the trail at the tip of mu fades out
SHOW_TIME_TRACKER = False # Not functional at the moment
ZOOM = 0.75 # Changes how close the camera is to the ball; good values for proper rolling are around 0.3, for oscillation around 0.75
ANIMATION_END = 15  # How many seconds of simulation time should be animated (e.g. if you simulated 20s, you may only want to show the interesting first 15 seconds)
RUN_TIME = 30 # Actual runtime of the final animation; if the ball moves very fast, factors of 20 or 30 * animation time are good; for small oscillations, *3 or *4 seems enough; for testing, use 1s
FRAME_TIME = 0.75 # Centeres the camera at the mean of the com's of all balls at time FRAME_TIME*(T_END-T_START)
POLAR = 65 # Polar angle of the camera in degree 
AZIMUTHAL = -45 # Azimuthal angle of the camera in degree
# ----------------------------
# Interpolators
# Do not change
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

# ----------------------------
# Animation
# If possible, try to change the global parameters instead for consistency
# ----------------------------

class Anim(ThreeDScene):
    def construct(self):
        time = ValueTracker(EXPORT_TIME if EXPORT_STILL else t[0])

        # frame_t determines the reference frame. It centeres the camera at the mean of all COM's after some portion of total time 
        frame_t = EXPORT_TIME if EXPORT_STILL else t[int(FRAME_TIME * (len(t) - 1))]
        all_points = np.array([state_at(frame_t, d)[0] for d in datasets])
        scene_center = SCALE * all_points.mean(axis=0)

        def world_pos(com):
            return SCALE * com - scene_center

        self.set_camera_orientation(
            phi=POLAR * DEGREES,
            theta=AZIMUTHAL * DEGREES,
            zoom=ZOOM,
        )

        balls = VGroup()
        labels = VGroup()
        trails = VGroup()
        mu_arrows = VGroup()
        mu_trails = VGroup()

        surface = Surface(
            lambda u, v: np.array([u, v, v * SLOPE]) * SCALE - scene_center,
            u_range=[SLOPE_X_MIN, SLOPE_X_MAX],
            v_range=[SLOPE_Y_MIN, SLOPE_Y_MAX],
            resolution=15,
            fill_opacity=SURFACE_OP,
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
                fill_opacity=SPHERE_OP,
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
            self.play(time.animate.set_value(ANIMATION_END), run_time=RUN_TIME, rate_func=linear)

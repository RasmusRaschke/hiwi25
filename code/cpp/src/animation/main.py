from manim import *
import numpy as np

phi = 0.175

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
            v_range=[-2, +2],
            u_range=[-2,+2]
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
        axes = ThreeDAxes()
        self.add(axes, plane, sphere)
        self.wait(10)
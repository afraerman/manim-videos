from manim import *
from manim.opengl import *
from scipy.integrate import odeint
from scipy.integrate import solve_ivp


def lorenz_system(t, state, sigma=10, rho=28, beta=8/3):
	x, y, z = state
	dxdt = sigma * (y - x)
	dydt = x * (rho - z) - y
	dzdt = x * y - beta * z

	return [dxdt, dydt, dzdt]

def ode_solution_points(function, state0, time, dt=0.01):
	solution = solve_ivp(
		function,
		t_span=(0, time),
		y0=state0,
		t_eval=np.arange(0, time, dt)
		)
	return solution.y.T


class LorenzAttractor(ThreeDScene):
	def construct(self):
		self.set_camera_orientation(phi=2*PI/5, theta=PI/5)
		axes = ThreeDAxes(
			x_range=(-50, 50, 5),
			y_range=(-50, 50, 5),
			z_range=(-0, 50, 5)
			)
		axes.center()
		self.add(axes)

		# Add the equations
		equations = Tex(
			R"""
			\begin{align}
            	\frac{\mathrm{d} x}{\mathrm{~d} t} & = \sigma(y-x)  \nonumber \\
            	\frac{\mathrm{d} y}{\mathrm{~d} t} & = x(\rho-z)-y \nonumber\\
            	\frac{\mathrm{d} z}{\mathrm{~d} t} & = x y-\beta z. \nonumber
		    \end{align}
			"""
		)
		xs = [1, 10, 18, 32]
		ys = [8, 13, 25, 33]
		zs = [22, 27, 36]
		for x in xs: equations[0][x].set_color(RED)
		for y in ys: equations[0][y].set_color(GREEN)
		for z in zs: equations[0][z].set_color(BLUE)

		equations.fix_in_frame()
		# equations.set_backstroke()
		equations.to_corner(UL)
		self.play(Write(equations))


		# Display Lorenz solutions
		epsilon = 1e-5
		number_of_lines = 10
		evolution_time = 30
		states = [
			[10, 10, 10+n*epsilon]
			for n in range(number_of_lines)
		]
		colors = color_gradient([BLUE_E, BLUE_A], number_of_lines)

		curves = OpenGLVGroup()

		for state, color in zip(states, colors):
			points = ode_solution_points(lorenz_system, state, evolution_time)

			curve = OpenGLVMobject()
			curve.set_points_as_corners(
				axes.c2p(points)
				)
			curve.set_stroke(color, 2)

			curves.add(curve)

		dots = Group(*(Dot3D(color=color) for color in colors))
		
		def update_dots(dots):
			for dot, curve in zip(dots, curves):
				dot.move_to(curve.get_end())

		dots.add_updater(update_dots)

		def update_curves(curves):
			number_of_points = 30 * evolution_time
			delta = 1 / (2/3 * number_of_points)
			curves.set_opacity(max(0, curves.get_opacity() - delta))

		curves.add_updater(update_curves)

		self.add(dots)
		self.add(curves)
		self.add(equations)
		self.play(
			*(
				Create(curve, rate_func=linear)
				for curve in curves
			),
			self.camera.animate.set_euler_angles(phi=2*PI/5, theta=PI),
			run_time=evolution_time,
		)
		self.interactive_embed()

class Hello(ThreeDScene):
	def construct(self):
		circle0 = Circle()
		circle = Circle()
		square = Square()
		circle.to_edge(LEFT)
		square.to_edge(UP)
		
		text = Text("Hello World", font_size=32)
		text.to_edge(UP)
		
		self.add(circle)
		self.add(square)
		
		self.play(
			Write(text, run_time=3),
			Transform(text[0], circle0),
			run_time=5,
			rate_func=linear
			
		)
		animations = [
			FadeOut(circle),
			FadeOut(square, shift=DOWN)
		]
		
		self.play(AnimationGroup(animations[0], lag_ratio=0.5))

		# self.interactive_embed()


import os
from numpy import std
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

plt.style.use("Solarize_Light2")


class Plotter:
    def __init__(
        self,
        n_states,
        title="Live Plot",
        xlabel="X-axis",
        ylabel="Y-axis",
        style="b-",
        save_path="animation.mp4",
    ):
        """Initializes the plotting class"""
        self.n_states = n_states
        self.x_data = []
        self.y_data = [[] for _ in range(self.n_states)]
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.style = style
        self.save_path = save_path

        self.fig, self.ax = plt.subplots()

        self.lines = []
        for i in range(self.n_states):
            """Create an empty plot for each state"""
            if self.title.__contains__("state"):
                self.lines.append(self.ax.plot([], [], label=f"x{i}"))
            else:
                """And for each input"""
                self.lines.append(self.ax.plot([], [], label=f"u{i}"))

        # Plot labels/details for each plot window

        self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.grid("on", color="black")
        self.ax.legend()

        self.writer = None

    def start_recording(self, fps=10):
        """Initialize plot recording."""
        self.writer = FFMpegWriter(
            fps=fps, metadata={"artist": "LivePlotter"}, bitrate=1800
        )
        self.writer.setup(self.fig, os.getcwd() + "/" + self.save_path, dpi=100)

    def update(self, x, y):
        """Update the plot with new x, y data points"""
        self.x_data.append(x)
        for i in range(self.n_states):
            self.y_data[i].append(y[i])

        state_bounds, input_bounds = None, None
        for i in range(self.n_states):
            if self.title.__contains__("state"):
                self.lines[i][0].set_data(self.x_data, self.y_data[i])
            else:
                self.lines[i][0].set_data(self.x_data, self.y_data[i])

        # Rescale plot y limits for each new y-appended value
        self.ax.relim()
        self.ax.autoscale_view()

        # Redraw
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        plt.pause(0.001)

        # Grab the current frame
        if self.writer is not None:
            self.writer.grab_frame()

    def stop_recording(self):
        """Finish and save"""
        if self.writer is not None:
            self.writer.finish()
            self.writer = None
            self.close()
            self.fig.savefig(f"{self.title}.png")

    def close(self):
        """Closes the plot."""
        plt.ioff()
        plt.close(self.fig)

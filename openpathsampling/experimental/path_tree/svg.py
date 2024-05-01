import xml.etree.ElementTree as ET
from .options import canonicalize_mover

def _stringify_values(dct):
    # also known as "XML is stupid"
    return {k: str(v) for k, v in dct.items()}

class SVGRendering:
    def __init__(self, hscale=10, vscale=10):
        self.hscale = hscale
        self.vscale = vscale
        self.trajectories = []
        self.connectors = []
        self.min_x = float("inf")
        self.max_x = float("-inf")

    @staticmethod
    def step_basics(step, options):
        mover = canonicalize_mover(step.mover)
        mover_options = options.movers[mover]
        color = mover_options.color
        plot_segments = mover_options.get_left_right(step)
        return mover, color, plot_segments

    def draw_trajectory(self, row, step, options):
        mover, color, plot_segments = self.step_basics(step, options)
        for left, right in plot_segments:
            if left < self.min_x:
                self.min_x = left
            if right > self.max_x:
                self.max_x = right

            attrib = {
                "width": (right - left) * self.hscale,
                "height": self.vscale // 2,
                "x": left * self.hscale,
                "y": row * self.vscale,
                "fill": color
            }
            elem = ET.Element("rect", attrib=_stringify_values(attrib))
            self.trajectories.append(elem)

    def draw_connector(self, x, bottom, top, step, options):
        attrib = {
            "x1": x * self.hscale,
            "y1": (bottom + 0.25) * self.vscale,
            "x2": x * self.hscale,
            "y2": (top + 0.5) * self.vscale,
            "style": "stroke:black",
        }
        elem = ET.Element("line", attrib=_stringify_values(attrib))
        self.connectors.append(elem)

    def build_svg(self):
        attrib = {
            "height": self.vscale * (len(self.trajectories) + 1),
            "width": self.hscale * (self.max_x - self.min_x),
            "xmlns": "http://www.w3.org/2000/svg"
        }
        root = ET.Element("svg", attrib=_stringify_values(attrib))
        for traj in self.trajectories:
            root.append(traj)
        for cnx in self.connectors:
            root.append(cnx)
        return root

    def draw(self):
        root = self.build_svg()
        return ET.tostring(root)

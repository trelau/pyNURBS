from pynurbs.geometry.creator import CreateGeom
from pynurbs.geometry.intersector import IntersectCurveCurve
from pynurbs.geometry.point import Point
from pynurbs.graphics.viewer import Viewer

gui = Viewer()

p0 = Point((0, 0, 0))
p1 = Point((2, 2, 0))
p2 = Point((4, 3, 0))
p3 = Point((8, -2, 0))

c1 = CreateGeom.interpolate_points([p0, p1, p2, p3])

gui.add_items(p0, p1, p2, p3, c1)
gui.show()

p0 = Point((5, -5, 0))
p1 = Point((10, 8, 0))

c2 = CreateGeom.interpolate_points([p0, p1])

gui.add_items(p0, p1, c2, c1)
gui.show()

cci = IntersectCurveCurve(c1, c2)

gui.add_items(c1, c2, *cci.points)
gui.show()

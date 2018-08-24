from pynurbs.geometry.creator import CreateGeom
from pynurbs.geometry.intersector import IntersectSurfaceSurface
from pynurbs.geometry.point import Point
from pynurbs.graphics.viewer import Viewer

gui = Viewer()

p0 = Point((0, 0, 0))
p1 = Point((2, 2, 0))
p2 = Point((4, 3, 0))
p3 = Point((8, -2, 0))
c1 = CreateGeom.interpolate_points([p0, p1, p2, p3])

p0 = Point((0, 0, 10))
p1 = Point((2, 2, 10))
p2 = Point((4, 3, 10))
p3 = Point((8, -2, 10))
c2 = CreateGeom.interpolate_points([p0, p1, p2, p3])

s1 = CreateGeom.interpolate_curves([c1, c2])

p0 = Point((0, 5, 5))
p1 = Point((2, 3, 5))
p2 = Point((4, 2, 5))
p3 = Point((8, 5, 5))
c3 = CreateGeom.interpolate_points([p0, p1, p2, p3])

p0 = Point((0, 5, 8))
p1 = Point((2, 3, 8))
p2 = Point((4, 2, 8))
p3 = Point((8, 5, 8))
c4 = CreateGeom.interpolate_points([p0, p1, p2, p3])

s2 = CreateGeom.interpolate_curves([c3, c4])

ssi = IntersectSurfaceSurface(s1, s2)

# Need Python implementation of surface tessellation to visualize...
gui.add_items(p0, p1, p2, p3, c1, c2, c3, c4, *ssi.icurves)
gui.show()

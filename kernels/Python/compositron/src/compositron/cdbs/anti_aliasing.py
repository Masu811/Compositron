from dataclasses import dataclass
from pathlib import Path

import ctypes
import numpy as np
import numpy.typing as npt

from ..constants import PI


agg = ctypes.CDLL(Path(__file__).parent / "../../../../../shared/agg/libagg.so")

agg.agg_aa.restype = None
agg.agg_aa.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.uint8, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    ctypes.c_int
]



@dataclass
class Vertex:
    x: float
    y: float


@dataclass
class Polygon:
    vertices: list[Vertex]

    def to_vertices(self) -> list[Vertex]:
        x0 = sum(vertex.x for vertex in self.vertices) / len(self.vertices)
        y0 = sum(vertex.y for vertex in self.vertices) / len(self.vertices)

        angles = [np.atan2(v.y - y0, v.x - x0) for v in self.vertices]

        return [self.vertices[i] for i in np.argsort(angles)]


@dataclass
class Rectangle:
    i_min: float
    i_max: float
    j_min: float
    j_max: float

    def to_vertices(self) -> list[Vertex]:
        return [
            Vertex(self.j_min, self.i_min),
            Vertex(self.j_min, self.i_max),
            Vertex(self.j_max, self.i_max),
            Vertex(self.j_max, self.i_min),
        ]


@dataclass
class Parallelogram:
    center_i: float
    center_j: float
    slope_1: float
    slope_2: float
    side_length_1: float
    side_length_2: float

    def to_vertices(self) -> list[Vertex]:
        half_width_1 = self.side_length_1 / 2
        half_width_2 = self.side_length_2 / 2

        w = np.array([1., self.slope_1])
        h = np.array([1., self.slope_2])

        w *= half_width_1 / np.linalg.norm(w)
        h *= half_width_2 / np.linalg.norm(h)

        x0 = np.array([self.center_j, self.center_i])

        v1 = x0 + (w + h)
        v2 = x0 + (w - h)
        v3 = x0 + (-w - h)
        v4 = x0 + (-w + h)

        return [
            Vertex(v1[0], v1[1]), Vertex(v2[0], v2[1]),
            Vertex(v3[0], v3[1]), Vertex(v4[0], v4[1]),
        ]


@dataclass
class Ellipse:
    center_i: float
    center_j: float
    radius_i: float
    radius_j: float
    phi: float

    def to_vertices(self) -> list[Vertex]:
        sin_phi = np.sin(self.phi)
        cos_phi = np.cos(self.phi)

        n_vertices = 128

        vertices = []

        for k in range(n_vertices):
            t = 2 * PI * k / n_vertices
            s = np.sin(t)
            c = np.cos(t)

            vertices.append(
                Vertex(
                    self.center_j + self.radius_i * c * sin_phi
                        + self.radius_j * s * cos_phi,
                    self.center_i + self.radius_i * c * cos_phi
                        - self.radius_j * s * sin_phi,
                )
            )

        return vertices


ConvexToVerticesCounterClockwise = Polygon | Rectangle | Parallelogram | Ellipse


def get_bounding_box(vertices: list[Vertex]) -> tuple[int, int, int, int]:
    left_col = int(round(min(v.x for v in vertices)))
    right_col = int(round(max(v.x for v in vertices)))
    lower_row = int(round(max(v.y for v in vertices)))
    upper_row = int(round(min(v.y for v in vertices)))

    return (left_col, right_col, lower_row, upper_row)


def agg_aa(
    nrows: int, ncols: int, polygon: ConvexToVerticesCounterClockwise
) -> npt.NDArray[np.uint8]:
    if nrows < 3 or ncols < 3:
        raise ValueError("Input data is too small in size")

    vertices = polygon.to_vertices()

    if len(vertices) < 3:
        raise ValueError("Insufficient number of vertices")

    weights = np.zeros((nrows, ncols), dtype=np.uint8)

    # Transform our coords (pixel centers at integer coords)
    # to agg's coords (pixel boundaries at integer coords)
    vertices = [Vertex(v.x + 0.5, v.y + 0.5) for v in vertices]

    if any(v.x < 0 or v.x > ncols or v.y < 0 or v.y > nrows for v in vertices):
        raise ValueError("Polygon vertices are out of matrix bounds")

    vx = np.asarray([v.x for v in vertices])
    vy = np.asarray([v.y for v in vertices])

    agg.agg_aa(weights, nrows, ncols, vx, vy, len(vertices))

    return weights


def compute_intersection(a: Vertex, b: Vertex, v1: Vertex, v2: Vertex) -> Vertex:
    m1 = (a.y - b.y) / (a.x - b.x)
    t1 = a.y - m1 * a.x

    m2 = (v1.y - v2.y) / (v1.x - v2.x)
    t2 = v1.y - m2 * v1.x

    x0 = (t2 - t1) / (m1 - m2)
    y0 = m1 * x0 + t1

    return Vertex(x0, y0)


def inside(v: Vertex, v1: Vertex, v2: Vertex) -> bool:
    a = np.array([v1.x - v2.x, v1.y - v2.y])
    b = np.array([v1.x - v.x, v1.y - v.y])

    return a[0] * b[1] - a[1] * b[0] < 0


def intersect_convex_polygons(
    a: ConvexToVerticesCounterClockwise,
    b: ConvexToVerticesCounterClockwise,
) -> Polygon:
    va = a.to_vertices()
    vb = b.to_vertices()

    if len(va) < 3 or len(vb) < 3:
        raise ValueError("Insufficient number of vertices")

    la = len(va)

    output = vb

    for i in range(la):
        edge_start = va[i]
        edge_end = va[(i + 1) % la]

        input = output
        output = []

        for i in range(len(input)):
            curr = input[i]
            prev = input[(i + len(input) - 1) % len(input)]

            curr_inside = inside(curr, edge_start, edge_end)
            prev_inside = inside(prev, edge_start, edge_end)

            if curr_inside:
                if not prev_inside:
                    output.append(compute_intersection(
                        prev, curr, edge_start, edge_end
                    ))
                output.append(curr);
            elif prev_inside:
                output.append(compute_intersection(
                    prev, curr, edge_start, edge_end
                ))

    return Polygon(output)

"""
R-Visibility polygons as python package. Also computes atomic visibility regions as described in Couto et al.
"""
# ._cgal_bindings will only exist after compilation.
from ._cgal_bindings import (
    FieldNumber,
    Point,
    Polygon,
    PolygonWithHoles,
    VisibilityPolygonCalculator,
    repair,
    AVP_Arrangement,
    Arrangement,
)

__all__ = [
    "FieldNumber",
    "Point",
    "Polygon",
    "PolygonWithHoles",
    "VisibilityPolygonCalculator",
    "repair",
    "AVP_Arrangement",
    "Arrangement",
]

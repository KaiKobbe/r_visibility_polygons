# r_visibility_polygons
This package provides functionality to compute r-visibility polygons for vertex guards in orthogonal polygons.

Consider an orthogonal polygon $\mathcal{P}$ with vertices $V(\mathcal{P})$.
The r-visibility polygon $Vis(g)$ of a guard $g \in V(\mathcal{P})$ is the set of all points $p \in \mathcal{P}$ such that the axis-parallel rectangle defined by $g$ and $p$ is fully contained in $\mathcal{P}$.

This model is a common way to define visibility—e.g., in art gallery problems—as it is more robust than defining visibility via an unobstructed line segment.
Moreover, r-visibility leads to one of the rare cases where the [art gallery problem can be solved in polynomial time](https://www.worldscientific.com/doi/10.1142/S0218195907002264?srsltid=AfmBOoqt3RYH032ohWCjp-ofxwXvRfHrhsqg6SboC2H6rDCO3cLojF2v).

This package has been developed as part of our [solver for the dispersive art gallery problem](https://github.com/KaiKobbe/dispersive_agp_solver).
Currently, the package does not include extensive sanity checks and offers only limited functionality.
Further features are planned for future versions.

## Installation
Since the code is writen in C++ and integrated in Python via PyBind11, a modern C++ compiler is necessary to use this package.
Installing should be easy to do via the following command:

```bash
pip install --verbose git+https://github.com/KaiKobbe/r_visibility_polygons
```

Note that the installation might take some time as CGAL will be locally installed.

## License
While this code is licensed under the MIT license, it has a dependency on
[CGAL](https://www.cgal.org/), which is licensed under GPL. This may have some
implications for commercial use.
/**
 * Specifies abstract classes and interfaces that define the nature of the physical space in which
 * a simulation is performed.  Implementations are given in the packages {@link etomica.space3d etomica.space3d},
 * {@link etomica.space2d etomica.space2d}, and {@link etomica.space1d etomica.space1d}, for 3-, 2-, and 1-dimensional
 * spaces, respectively.  The key elements are:
 * <ul>
 * <li>Space, which serves primarily as a factory for the other classes.  Subclasses of Space make vectors, etc.
 * appropriate to the space it represents.
 * <li>Vector, which for a D-dimensional space is a set of D fields of type double, with methods used the
 * perform vector operations on them.
 * <li>Tensor, which defines a 2nd-rank tensor in the D-dimensional space.  It hold D<sup>2</sup> fields of type double, and
 * defines methods used to perform tensor operations on them.
 * <li>Orientation, which provides fields and operations define and manipulate the orientation of a rigid body.
 * <li>Boundary, which is held by a Box instance to define the size of the box and to specify what
 * happens as atoms cross the boundary of the space (e.g., periodic boundaries).
 * </ul>
 * All of the classes above (except perhaps Boundary) are implemented in the space subpackages in a manner
 * specific to the space being defined.
 *
 */
package etomica.space;
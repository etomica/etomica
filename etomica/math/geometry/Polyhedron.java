package etomica.math.geometry;
import etomica.*;
import etomica.space3d.Space3D;

/**
 * Representation of a mathematical polygon, a polytope in 3-D.
 */
public abstract class Polyhedron extends Polytope {
    
    public Polyhedron() {
        super(Space3D.INSTANCE);
    }
     
 }

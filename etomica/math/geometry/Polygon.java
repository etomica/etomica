package etomica.math.geometry;
import etomica.*;
import etomica.space2d.Space2D;

/**
 * Representation of a mathematical polygon, a polytope in 2-D.
 */
public abstract class Polygon extends Polytope {
    
    public Polygon() {
        super(Space2D.INSTANCE);
    }
     
 }

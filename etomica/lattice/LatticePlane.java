package etomica.lattice;

import etomica.math.geometry.Plane;
import etomica.*;

/**
 * Class describing a plane through a lattice.
 */
 
 /* History
  * 09/07/02 (DAK) new
  */
  
public class LatticePlane implements AtomFilter {
    
    private final Plane plane;
    
    public LatticePlane(int h, int k, int l) {
        plane = Plane.newInterceptForm(h, k, l);
    }
    
    public boolean accept(Atom a) {
        return plane.isPositiveSide((Space3D.Vector)a.coord.position());
    }
    
    public Plane getPlane() {return plane;}
    
}//end of LatticePlane
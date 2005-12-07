/*
 * Created on Dec 6, 2004
 *
 */
package etomica.models.hexane;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.space.Space;

/**
 * @author nancycribbin
 * Implements the triclinic crystal structure derived from Dr. Monson's data.
 */

public class PrimitiveHexane extends PrimitiveTriclinic {
    
    public PrimitiveHexane(Space space) {
        //In units of sigma
        super(space, 2.86842545, 1.338313769, 1.290909603, 0.779541414, 1.109209538,
                0.768992016);
        
        //In units of Angstroms
//        super(space, 28.68425449956680, 13.38313769270750, 12.90909602969220, 
//                0.779541414, 1.109209538, 0.768992016);
    }
}

//double a, double b, double c, double alpha, double beta, double gamma)

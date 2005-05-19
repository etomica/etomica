/*
 * Created on Dec 6, 2004
 *
 */
package etomica.models.hexane;
import etomica.Space;
import etomica.lattice.crystal.PrimitiveTriclinic;

/**
 * @author nancycribbin
 * Implements the triclinic crystal structure derived from Dr. Monson's data.
 */

public class PrimitiveHexane extends PrimitiveTriclinic {
    
    public PrimitiveHexane(Space space) {
        super(space, 2.86842545, 1.338313769, 1.290909603, 0.779541414, 1.109209538,
                0.768992016);
    }
}

//double a, double b, double c, double alpha, double beta, double gamma)

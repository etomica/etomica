package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.potential.IPotentialAtomicMultibody;

/**
 * Non-additive Mayer function class for "spherical" potentials 
 *
 * @author Andrew Schultz
 */
public class MayerFunctionSphericalThreeBody extends MayerFunctionThreeBody {

    protected final IPotentialAtomicMultibody p3;
    
    public MayerFunctionSphericalThreeBody(IPotentialAtomicMultibody p3) {
        this.p3 = p3;
    }
    
    protected double energy(IMoleculeList molecules, double[] r2) {
        return p3.energy(r2);
    }

    public void setBox(IBox box) {
        p3.setBox(box);
        super.setBox(box);
    }
}

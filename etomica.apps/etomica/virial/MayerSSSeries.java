package etomica.virial;

import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.space.ISpace;

/**
 * @author kofke
 *
 * Hard-sphere Mayer function.  -1 if r < sigma; 0 otherwise
 */
public class MayerSSSeries extends MayerFunctionSpherical {

    private static final long serialVersionUID = 1L;
    protected int exp6;

    /**
     * Constructor for MayerHardSphere.
     */
    public MayerSSSeries(ISpace _space, int exp) {
        super(_space);
        setExp(exp);
    }

    /**
     * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair)
     */
    public double f(double r2, double beta) {
        // beta = 1, epsilon = 1, sigma = 1
        if (r2 == 0) {
            return 0;
        }
        double s6 = 1.0/(r2*r2*r2);
        double sexp = s6;
        for (int i=1; i<exp6; i++) {
            sexp *= s6;
        }
        return Math.exp(-(4.0*s6*s6))*sexp;
    }

    /**
     * Returns the HS diameter.
     * @return double
     */
    public int getExp() {
        return exp6*6;
    }

    /**
     * Sets the HS diameter.
     * @param sigma The sigma to set
     */
    public void setExp(int newExp) {
        exp6 = newExp/6;
    }

    public void setBox(IBox newBox) {
    }
    
    public IPotential getPotential() {
        return null;
    }
}

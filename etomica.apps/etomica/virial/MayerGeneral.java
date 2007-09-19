package etomica.virial;

import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.potential.IPotential;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerGeneral implements MayerFunction, java.io.Serializable {

    /**
     * Constructor Mayer function using given potential.
     */
    public MayerGeneral(IPotential potential) {
        this.potential = potential;
    }

    public double f(AtomSet pair, double beta) {
        double betaU = beta*potential.energy(pair);
        if (Math.abs(betaU) < 1.e-8) {
            // for small betaU, exp(-betaU)-1 ~= -betaU
            // for betaU < 1E-8, the approximation is value within machine precision
            // for betaU < 1E-15, exp(-betaU) is 1, so the approximation is more accurate
            //   than simply doing the math.
            return -betaU;
        }
        return Math.exp(-betaU) - 1.0;
    }

    public IPotential getPotential() {
        return potential;
    }

    public void setBox(Box newBox) {
        potential.setBox(newBox);
    }

    private final IPotential potential;
}

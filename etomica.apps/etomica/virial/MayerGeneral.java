package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.api.IPotentialMolecular;

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
    public MayerGeneral(IPotentialMolecular potential) {
        this.potential = potential;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        double x = -beta*potential.energy(pair);
        if (Math.abs(x) < 0.01) {
            return x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        return Math.exp(x) - 1;
    }

    public IPotential getPotential() {
        return potential;
    }

    public void setBox(IBox newBox) {
        potential.setBox(newBox);
    }

    private final IPotentialMolecular potential;
}

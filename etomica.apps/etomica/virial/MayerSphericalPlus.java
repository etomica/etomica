package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.potential.Potential2Spherical;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerSphericalPlus extends MayerGeneralSpherical {

    /**
     * Constructor Mayer function using given potential.
     */
    public MayerSphericalPlus(Potential2Spherical potential, Potential2Spherical potentialPlus) {
        super(potential);
        this.potentialPlus = potentialPlus;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        double f = super.f(pair, r2, beta);
        f += 2-2*Math.exp(-potentialPlus.u(r2)*beta);
        return f;
    }

    public IPotential getPotentialPlus() {
        return potentialPlus;
    }

    public void setBox(IBox newBox) {
        super.setBox(newBox);
        potentialPlus.setBox(newBox);
    }

    protected final Potential2Spherical potentialPlus;
}

// This class collects all the quantities that might be of interest
// when potential computes the interactions between a pair of Elements.
// All quantities can be computed at once almost as easily as any one
// of them, so by returning an instance of this class all quantities
// become available without having to recompute them separately

package simulate;

public class PairInteraction {

    public double[] force, rij;
    public double energy, virial, rSquared, dfdr;

    public PairInteraction() {
        force = new double[Space.D];
        rij = new double[Space.D];
        Space.uEa1(force,0.0);
        Space.uEa1(rij,0.0);
        rSquared = Double.MAX_VALUE;
        energy = virial = 0.0;
    }

}


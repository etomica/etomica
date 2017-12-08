package etomica.potential;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.numerical.BesselFunction;

public class PotentialCalculationFSum implements PotentialCalculationMolecular {
    protected final Vector ei, ej;
    protected final double Q, mu, J, bt; //TODO should I add final here
    protected Vector dr;
    protected double FSum = 0;
    protected DipoleSource dipoleSource;


    public PotentialCalculationFSum(Space space, double dipoleMagnitude, double interactionS, double temperature) {
        dr = space.makeVector();
        ei = space.makeVector();
        ej = space.makeVector();

        J = interactionS;
        mu = dipoleMagnitude;
        bt = 1 / temperature;
        Q = bt * bt * mu * mu * (1 + BesselFunction.I(1, J * bt) / BesselFunction.I(0, J * bt));
    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {//TODO ??? should I do similar stuff here????
            return;
        }
        //get total force here!!!
        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;

        Vector[][] t = potentialSecondDerivative.gradientAndTorque(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);

        //TODO the acos returns 0 to Pi but t1 is form 0 to 2Pi
        double x1 = atom1.getOrientation().getDirection().getX(0);//cost1
        double y1 = atom1.getOrientation().getDirection().getX(1);//sint1
        double x2 = atom2.getOrientation().getDirection().getX(0);//cost2
        double y2 = atom2.getOrientation().getDirection().getX(1);//sint2
        if (x1 > 1) {
            x1 = 1;
        } else if (x1 < -1) {
            x1 = -1;
        }
        if (y1 > 1) {
            y1 = 1;
        } else if (y1 < -1) {
            y1 = -1;
        }
        if (x2 > 1) {
            x2 = 1;
        } else if (x2 < -1) {
            x2 = -1;
        }
        if (y2 > 1) {
            y2 = 1;
        } else if (y2 < -1) {
            y2 = -1;
        }
        double t1 = Math.acos(x1);
        double t2 = Math.acos(x2);
        double bt2 = bt * bt;
        double bt3 = bt2 * bt;
        double mu2 = mu * mu;
        double mu3 = mu2 * mu;//TOOD check what you need to return, sum F * sum FC or other ways

        //Try to implement the constant part here
        FSum += -2 * Q + J * bt * bt2 * (t1 - t2) * mu2 * (Math.sin(2 * t2) - Math.sin(2 * t1));

    }

    public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
        if (!(potential instanceof IPotentialMolecularSecondDerivative)) {
            return;
        }


    }

    public void setDipoleSource(DipoleSource newDipoleSource) {
        dipoleSource = newDipoleSource;
    }

    public void zeroSum() {
        FSum = 0.0;
    }

    /**
     * Returns the current value of the energy sum.
     */
    public double getSum() {
        return FSum;
    }


}

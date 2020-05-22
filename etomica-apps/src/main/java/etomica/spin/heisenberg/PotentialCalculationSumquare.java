package etomica.spin.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Anx expressions for computing v_E and v_EE in the mapping.
 *
 * @author Weisong Lin
 */

public class PotentialCalculationSumquare implements PotentialCalculation {
    protected Vector ei, ej;
    protected double sum, mu;
    protected Vector dr;

    public PotentialCalculationSumquare(Space space, double dipoleMagnitude) {
        ei = space.makeVector();
        ej = space.makeVector();
        dr = space.makeVector();
        mu = dipoleMagnitude;

    }


    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());

//        System.out.println("conventional pairs:("+ atom1+" , "+ atom2+")");

        dr.E(ei);
        dr.PE(ej);

        sum += mu * mu * dr.squared();
    }

    public void zeroSum() {
        sum = 0;
    }

    public double getSum() {
        return sum;
    }


}

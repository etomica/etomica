package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.util.numerical.BesselFunction;

public class PotentialCalculationPhiSumHeisenberg implements PotentialCalculationMolecular {
    protected Vector ei, ej;
    protected Vector dr;
    protected double secondDerivativeSum = 0;
    protected DipoleSource dipoleSource;
    protected double bt, J;


    public PotentialCalculationPhiSumHeisenberg(Space space, double interactionS, double beta) {
        dr = space.makeVector();
        ei = space.makeVector();
        ej = space.makeVector();
        J = interactionS;
        bt = beta;

    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }

        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());
//try to input all the function here,please!
        double bJ = bt * J;
        double cos = ei.dot(ej);
        double exp2bJ = Math.exp(-2 * bJ * cos);
//        double bessel0 = BesselFunction.doCalc(true,0,bJ);
//        double bessel1 = BesselFunction.doCalc(true,1,bJ);
        double bessel = BesselFunction.doCalc(true, 0, bJ)
                + BesselFunction.doCalc(true, 1, bJ);

        secondDerivativeSum += t[1].component(0, 0) * exp2bJ * bessel * bessel;
        //t[1] is ii component! p11,note it's not times by bt yet, should times it back in meter!!
    }

    public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
        if (!(potential instanceof IPotentialMolecularSecondDerivative)) {
            return;
        }
    }


    public void zeroSum() {
        secondDerivativeSum = 0.0;
    }

    /**
     * Returns the current value of the energy sum.
     */
    public double getSum() {
        return secondDerivativeSum;
    }


}

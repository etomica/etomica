package etomica.potential;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCalculationPhiSumHeisenberg implements PotentialCalculationMolecular {
    protected Vector ei, ej;
    protected Vector dr;
    protected double secondDerivativeSum = 0;


    public PotentialCalculationPhiSumHeisenberg(Space space) {
        dr = space.makeVector();
        ei = space.makeVector();
        ej = space.makeVector();

    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }

        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.get(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.get(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());

        secondDerivativeSum += 2 * t[0].component(0, 0) * (ei.dot(ej));
        secondDerivativeSum += t[1].component(0, 0);
        secondDerivativeSum += t[2].component(0, 0);

//		double Cos = ei.dot(ej);
//		secondDerivativeSum += -2*Cos*Cos+2*Cos;
//		double diff = 2*t[0].component(0,0)*(ei.dot(ej))+t[1].component(0,0)+t[2].component(0,0)-1.5*( -2*Cos*Cos+2*Cos);
//		System.out.println("diff = " + diff);//check these two approach is the same or not
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

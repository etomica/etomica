package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomOriented;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.util.numerical.BesselFunction;

public class PotentialCalculationPhiSumHeisenberg implements PotentialCalculationMolecular {
    protected final IVectorMutable ei, ej;
    protected IVectorMutable dr;
    protected double secondDerivativeSum = 0;
    protected DipoleSource dipoleSource;


    public PotentialCalculationPhiSumHeisenberg(ISpace space) {
        dr = space.makeVector();
        ei = space.makeVector();
        ej = space.makeVector();

    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());

        //complicated way to get secondDerivativeSum, don't need the secondDerivative from p2Spin
        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);
        double Cos = ei.dot(ej);
        secondDerivativeSum += 2.0 * t[0].component(0, 0) * Cos;
        secondDerivativeSum += t[1].component(0, 0);
        secondDerivativeSum += t[2].component(0, 0);

//        System.out.println("ei in PhiSum = " + ei);
//        System.out.println(t[0].component(0, 0));
//        System.exit(2);

        //much easier way but need to time back coupling J in the meter and here!!!!!!!!
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

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
	 protected final IVectorMutable ei,ej;
	 protected IVectorMutable dr;
	 protected double secondDerivativeSum= 0;
	 protected DipoleSource dipoleSource;
	 
	 protected  double Q ,mu,J,bt;

	 
	public PotentialCalculationPhiSumHeisenberg(ISpace space) {
	    dr = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	    
	}

	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		if(!(potential instanceof IPotentialAtomicSecondDerivative)){
			return;
		}
		IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;

		Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);
		
		IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(0);
    	IAtomOriented atom2 = (IAtomOriented)atoms.getAtom(1);
    	
    	ei.E(atom1.getOrientation().getDirection());
    	ej.E(atom2.getOrientation().getDirection());
    	

		double s1 = ei.getX(1);
		double s2 = ej.getX(1);


		secondDerivativeSum += 2*t[0].component(0, 0)*s1*s2 + t[1].component(0, 0)*s1*s1
				+t[2].component(0, 0)*s2*s2;

//		System.out.println( t[0].component(0, 0) + " " + t[1].component(0, 0) + " " + t[2].component(0, 0));
//		System.out.println("secondDerivative = " + secondDerivativeSum);
//		System.exit(2);
	}

	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if(!(potential instanceof IPotentialMolecularSecondDerivative)){
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

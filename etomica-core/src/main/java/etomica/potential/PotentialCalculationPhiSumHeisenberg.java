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
	 

	public PotentialCalculationPhiSumHeisenberg(ISpace space) {
	    dr = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	    
	}

	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		if(!(potential instanceof IPotentialAtomicSecondDerivative)){
			return;
		}

		//don't need this for now
//		IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
//
//		Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);
		
		IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(0);
    	IAtomOriented atom2 = (IAtomOriented)atoms.getAtom(1);
    	ei.E(atom1.getOrientation().getDirection());
    	ej.E(atom2.getOrientation().getDirection());

		double c1 = ei.getX(0);//cost1
		double c2 = ej.getX(0);//cost2
		double s1 = ei.getX(1);//sint1
		double s2 = ej.getX(1);//sint2

		//ij phi_ij = -J*(1+c1*c2/s1/s2); J would be multiplied in meter
		//-J Sin[t1] Sin[t2] (Cos[t1] Cos[t2] + Sin[t1] Sin[t2])
//		secondDerivativeSum += -2.0*(s1*s2+c1*c2)*s1*s2;

		//ii phi_ii = J*s2^3/s1; J would be multiplied in meter
		//J Sin[t1] Sin[t2]
		//jj phi_jj = J*s1^3/s2; J would be multiplied in meter
		//J Sin[t1] Sin[t2]
		//ii and jj is the same
//		secondDerivativeSum += 2.0*s1*s2;

		//or you could combine ij ii and jj
//		secondDerivativeSum += -2.0*(s1*s2+c1*c2-1)*s1*s2;

		double Cos = ei.dot(ej);
		secondDerivativeSum += -2*Cos*Cos+2*Cos;
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

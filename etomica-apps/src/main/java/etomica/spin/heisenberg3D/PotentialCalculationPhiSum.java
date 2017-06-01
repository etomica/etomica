package etomica.spin.heisenberg3D;

import etomica.api.*;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.IPotentialMolecularSecondDerivative;
import etomica.potential.PotentialCalculationMolecular;
import etomica.space.ISpace;
import etomica.space.Tensor;

public class PotentialCalculationPhiSum implements PotentialCalculationMolecular {
	 protected IVectorMutable fieldE;
	 protected final IVectorMutable ei,ej;
	 protected IVectorMutable Ai;
	 protected IVectorMutable Aj;
	 protected IVectorMutable dr;
	 protected double secondDerivativeSum= 0;
	 protected DipoleSource dipoleSource;
	 protected final IVectorMutable [] a;
	 protected final Tensor iT;
	 
	public PotentialCalculationPhiSum(ISpace space) {
		fieldE = space.makeVector();
		Ai = space.makeVector();
	    Aj = space.makeVector();
	    dr = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	    a = new IVectorMutable[3];
		a[0] = space.makeVector();
		a[1] = space.makeVector();
		a[2] = space.makeVector();
		iT = space.makeTensor();
		double [] xD = {1,0,0};
		double [] yD = {0,1,0};
		double [] zD = {0,0,1};
		a[0].E(xD);
		a[1].E(yD);
		a[2].E(zD);
		iT.E(a);
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

		double c1 = ei.getX(0);//cost1
		double c2 = ej.getX(0);//cost2
		double s1 = ei.getX(1);//sint1
		double s2 = ej.getX(1);//sint2

		double traceij = t[0].trace();
		double traceii = t[1].trace();
		double tracejj = t[2].trace();


		t[0].transpose();
		t[0].TE(-1);
		t[1].transpose();
		t[1].TE(-1);
		t[2].transpose();
		t[2].TE(-1);


		t[0].PEa1Tt1(traceij, iT);
		t[1].PEa1Tt1(traceii, iT);
		t[2].PEa1Tt1(tracejj, iT);



		dr.E(ej);
		t[0].transform(dr);
		secondDerivativeSum += 2*ei.dot(dr);//ij

		dr.E(ei);
		t[1].transform(dr);
		secondDerivativeSum += ei.dot(dr);//ii

		dr.E(ej);
		t[2].transform(dr);
		secondDerivativeSum += ej.dot(dr);//jj

	}

	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if(!(potential instanceof IPotentialMolecularSecondDerivative)){
			return;
		}
		
		IPotentialMolecularSecondDerivative potentialSecondDerivative = (IPotentialMolecularSecondDerivative) potential;
		
		Tensor[] t = potentialSecondDerivative.secondDerivative(molecules);
		
		IMolecule molecule0 = molecules.getMolecule(0);
		IMolecule molecule1 = molecules.getMolecule(1);
		ei.E(dipoleSource.getDipole(molecule0));
		ej.E(dipoleSource.getDipole(molecule1));
		ei.normalize();
		ej.normalize();

		double traceij = t[0].trace();
		double traceii = t[1].trace();
		double tracejj = t[2].trace();
		
		
		t[0].transpose();
		t[0].TE(-1);
		t[1].transpose();
		t[1].TE(-1);
		t[2].transpose();
		t[2].TE(-1);

	
		t[0].PEa1Tt1(traceij, iT);
		t[1].PEa1Tt1(traceii, iT);
		t[2].PEa1Tt1(tracejj, iT);
		

		
		dr.E(ej);
		t[0].transform(dr);
		secondDerivativeSum += 2*ei.dot(dr);//ij

		
		
		dr.E(ei);
		t[1].transform(dr);
		
		secondDerivativeSum += ei.dot(dr);//ii
		
		dr.E(ej);
		t[2].transform(dr);
		
		secondDerivativeSum += ej.dot(dr);//jj
		

		
	}
	
	 public void setDipoleSource(DipoleSource newDipoleSource) {
	        dipoleSource = newDipoleSource;
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

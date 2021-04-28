/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.molecule.DipoleSourceMolecular;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCalculationPhiSum implements PotentialCalculationMolecular {
	 protected Vector fieldE;
	protected Vector ei, ej;
	 protected Vector Ai;
	 protected Vector Aj;
	 protected Vector dr;
	 protected double secondDerivativeSum= 0;
	 protected DipoleSourceMolecular dipoleSource;
	 protected final Vector[] a;
	 protected final Tensor iT;
	 
	public PotentialCalculationPhiSum(Space space) {
		fieldE = space.makeVector();
		Ai = space.makeVector();
	    Aj = space.makeVector();
	    dr = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	    a = new Vector[3];
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

	}

	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if(!(potential instanceof IPotentialMolecularSecondDerivative)){
			return;
		}
		
		IPotentialMolecularSecondDerivative potentialSecondDerivative = (IPotentialMolecularSecondDerivative) potential;
		
		Tensor[] t = potentialSecondDerivative.secondDerivative(molecules);
		
		IMolecule molecule0 = molecules.get(0);
		IMolecule molecule1 = molecules.get(1);
		ei.E(dipoleSource.getDipole(molecule0));
		ej.E(dipoleSource.getDipole(molecule1));
		ei.normalize();
		ej.normalize();
//		System.out.println("ei = " + ei);
//		System.out.println("ej = " + ej);
//		System.exit(2);
		
		
//		ei.normalize();
//		ej.normalize();
		
//		debug only  
//		IVectorMutable pos0 = atom1.getPosition();
//		IVectorMutable pos1 = atom0.getPosition();
//		Ai.E(pos1);
//		Ai.ME(pos0);
//		System.out.println("r = " + Ai);
//		System.out.println("ei = " + ei);
//		System.out.println("ej = " + ej);
		
		
		
//		System.out.println("ei = " + ei);
//		System.out.println("ej = " + ej);
//		System.out.println("t[0] = \n" + t[0]);
//		System.out.println("trace t[0] = " + t[0].trace());
		
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
		
//		System.out.println("-Transpose(t[0]) + trace(t[0]) = \n"  + t[0]);
//		System.exit(2);
		
		
		dr.E(ej);
		t[0].transform(dr);
		secondDerivativeSum += 2*ei.dot(dr);//ij
		
//		System.out.println("ij = " + 2*ei.dot(dr));
//		System.exit(2);
		
		
		
		dr.E(ei);
		t[1].transform(dr);
		
		secondDerivativeSum += ei.dot(dr);//ii
		
		dr.E(ej);
		t[2].transform(dr);
		
		secondDerivativeSum += ej.dot(dr);//jj
		

		
	}
	
	 public void setDipoleSource(DipoleSourceMolecular newDipoleSource) {
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

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
	 
	 protected  double Q ,mu,J,bt; //TODO shoud I add final?? 

	 
	public PotentialCalculationPhiSumHeisenberg(ISpace space, double dipoleMagnitude, double interactionS, double temperature) {
	    dr = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	    
		J = interactionS;
		mu = dipoleMagnitude;
		bt = 1/temperature;
		
    	Q = bt*bt*mu*mu*(1+BesselFunction.I(1, J*bt)/BesselFunction.I(0, J*bt));
//    	System.out.println("J="+ J+  " mu= "+ mu +" bt= "+ bt );
//    	System.out.println("Q= " + Q);
		
	}

	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		if(!(potential instanceof IPotentialAtomicSecondDerivative)){
			return;
		}
		IPotentialAtomicSecondDerivative potentialSeconDerivative = (IPotentialAtomicSecondDerivative) potential;
		
		
		Tensor[] t = potentialSeconDerivative.secondDerivative(atoms);
		
		IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(0);
    	IAtomOriented atom2 = (IAtomOriented)atoms.getAtom(1);
    	
//		ei.E(dipoleSource.getDipole(atom1.getParentGroup()));
//		ej.E(dipoleSource.getDipole(atom2.getParentGroup()));
    	ei.E(atom1.getOrientation().getDirection());
    	ej.E(atom2.getOrientation().getDirection());
    	
		

		if(ei.getX(0) > 1){
			ei.setX(0, 1);
		}
		if(ei.getX(0)<-1){
			ei.setX(0, -1);
		}
		if(ej.getX(0) > 1){
			ej.setX(0, 1);
		}
		if(ej.getX(0)<-1){
			ej.setX(0, -1);
		}
		
		double t1 = Math.acos(ei.getX(0));
		double t2 = Math.acos(ej.getX(0));
		
		double bt2=bt*bt;
		double bt3=bt*bt*bt;
		double mu2=mu*mu;
		double phiC = -2*Q+2*bt2*mu2+2*bt2*mu2*Math.cos((t1-t2))-J*bt3*(t1-t2)*mu2*Math.sin(2*t1)
				     +J*bt3*t1*mu2*Math.sin(2*t2)-J*bt3*t2*mu2*Math.sin(2*t2);
		
		
		
		//is this += or just =??????????? TODO also this part is always zero!!!!!!!!!!!!!
	
		secondDerivativeSum += t[0].component(0, 0)*phiC+t[1].component(0, 0)*phiC
				+ t[1].component(0, 0)*phiC+t[2].component(0, 0)*phiC;	
//		System.out.println("secondDerivative = " + secondDerivativeSum);
//		System.exit(2);
	}

	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if(!(potential instanceof IPotentialMolecularSecondDerivative)){
			return;
		}
		
//		IPotentialMolecularSecondDerivative potentialSeconDerivative = (IPotentialMolecularSecondDerivative) potential;
//		
//		Tensor[] t = potentialSeconDerivative.secondDerivative(molecules);
//		
//		IMolecule molecule0 = molecules.getMolecule(0);
//		IMolecule molecule1 = molecules.getMolecule(1);
//		
//		ei.E(dipoleSource.getDipole(molecule0));
//		ej.E(dipoleSource.getDipole(molecule1));
//		ei.normalize();
//		ej.normalize();
//		
//		System.out.println("ei = " + ei);
//		System.out.println("ej = " + ej);
//		System.exit(2);
//		
//		ei.normalize();
//		ej.normalize();
//		
//		if(ei.getX(0) > 1){
//			ei.setX(0, 1);
//		}
//		if(ei.getX(0)<-1){
//			ei.setX(0, -1);
//		}
//		if(ej.getX(0) > 1){
//			ej.setX(0, 1);
//		}
//		if(ej.getX(0)<-1){
//			ej.setX(0, -1);
//		}
//		
//		double t1 = Math.acos(ei.getX(0));//TODO anything I shoul be aware??? still need to handle abs>1 case
//		double t2 = Math.acos(ej.getX(0));
//		
//		double bt2=bt*bt;
//		double bt3=bt*bt*bt;
//		double mu2=mu*mu;
//		double phiC = -2*Q+2*bt2*mu2+2*bt2*mu2*Math.cos((t1-t2))-J*bt3*(t1-t2)*mu2*Math.sin(2*t1)
//				     +J*bt3*t1*mu2*Math.sin(2*t2)-J*bt3*t2*mu2*Math.sin(2*t2);
//		
//		
//		
////		System.out.println("J="+ J+  " mu= "+ mu +" bt= "+ bt + " phi= " + phiC);
//		secondDerivativeSum = t[0].component(0, 0)*phiC+t[1].component(0, 0)*phiC
//						 	+ t[1].component(0, 0)*phiC+t[2].component(0, 0)*phiC;
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

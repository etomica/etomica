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

public class PotentialCalculationFSum implements PotentialCalculationMolecular {
	 protected final IVectorMutable ei,ej;
	 protected IVectorMutable dr;
	 protected double FSum = 0;
	 protected DipoleSource dipoleSource;
	 
	 protected final double Q ,mu,J,bt; //TODO should I add final here

	 
	public PotentialCalculationFSum(ISpace space, double dipoleMagnitude, double interactionS, double temperature) {
	    dr = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	    
		J = interactionS;
		mu = dipoleMagnitude;
		bt = 1/temperature;
    	Q = bt*bt*mu*mu*(1+BesselFunction.I(1, J*bt)/BesselFunction.I(0, J*bt));
	}

	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		if(!(potential instanceof IPotentialAtomicSecondDerivative)){//TODO ??? should I do similar stuff here????
			return;
		}
		//get total force here!!!
		IPotentialAtomicSecondDerivative potentialSeconDerivative = (IPotentialAtomicSecondDerivative) potential;
		
		IVector[][] t = potentialSeconDerivative.gradientAndTorque(atoms);
		
		IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(0);
		IAtomOriented atom2 = (IAtomOriented)atoms.getAtom(1);
		
		//TODO the acos returns 0 to Pi but t1 is form 0 to 2Pi
		double x1 = atom1.getOrientation().getDirection().getX(0);//cost1
		double y1 = atom1.getOrientation().getDirection().getX(1);//sint1
		double x2 = atom2.getOrientation().getDirection().getX(0);//cost2
		double y2 = atom2.getOrientation().getDirection().getX(1);//sint2
		if(x1>1){x1=1;}else if(x1<-1){x1=-1;}
		if(y1>1){y1=1;}else if(y1<-1){y1=-1;}
		if(x2>1){x2=1;}else if(x2<-1){x2=-1;}
		if(y2>1){y2=1;}else if(y2<-1){y2=-1;}
		double t1 = Math.acos(x1);
		double t2 = Math.acos(x2);
		double bt2 = bt*bt;
		double bt3 = bt2*bt;
		double mu2 = mu*mu;
		double mu3 = mu2*mu;//TOOD check what you need to return, sum F * sum FC or other ways
		double f1C = 0.5*(-4*Q-2*Q*t1-2*Q*t2+2*bt2*mu2*t1+2*bt2*mu2*t2
				+J*bt3*(t1-t2)*mu2*Math.cos(2*t1)
				+2*bt2*(t1+t2)*mu2*Math.cos(t1-t2)-J*bt3*t1*mu2*Math.cos(2*t2)
				+J*bt3*t2*mu2*Math.cos(2*t2)-2*bt2*t1*mu2*Math.cos(t1+t2)
				+2*bt2*t2*mu2*Math.cos(t1+t2)-2*J*bt3*t1*mu2*Math.sin(2*t1)
				+2*J*bt3*t2*mu2*Math.sin(2*t1)+2*J*bt3*t1*mu2*Math.sin(2*t2)
				-2*J*bt3*t2*mu2*Math.sin(2*t2)+2*bt2*mu2*Math.sin(t1+t2));
		
		double f2C = 0.5*(  -4*Q-2*Q*t1-2*Q*t2+4*bt2*mu2 +2*bt2*mu2*t1+2*bt2*mu2*t2
				+J*bt3*(t1-t2)*mu2*Math.cos(2*t1)
				+2*bt2*(2+t1+t2)*mu2*Math.cos(t1-t2)-J*bt3*t1*mu2*Math.cos(2*t2)
				+J*bt3*t2*mu2*Math.cos(2*t2)+2*bt2*t1*mu2*Math.cos(t1+t2)
				-2*bt2*t2*mu2*Math.cos(t1+t2)-2*J*bt3*t1*mu2*Math.sin(2*t1)
				+2*J*bt3*t2*mu2*Math.sin(2*t1)+2*J*bt3*t1*mu2*Math.sin(2*t2)
				-2*J*bt3*t2*mu2*Math.sin(2*t2)+2*bt2*mu2*Math.sin(t1+t2) 	);

		double f1 = -J*Math.sin(t1-t2);//TODO valus not the same with torque in p2spin
		double f2 = -f1;
		System.out.println("=======================");
		System.out.println(" t1= "+t1);
		System.out.println("x1 = " + x1 + " y1 = " + y1);
		System.out.println("x2 = " + x2 + " y2 = " + y2);
		System.out.println("f1 = " + f1);
		System.out.println("f2 = " + f2);
		FSum += f1*f1C + f2*f2C;
//		System.exit(2);
		
		//Try to implement the constant part here
		FSum += -2*Q+J*bt*bt2*(t1-t2)*mu2*(Math.sin(2*t2)-Math.sin(2*t1));
		
	}

	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if(!(potential instanceof IPotentialMolecularSecondDerivative)){
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

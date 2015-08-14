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

public class PotentialCalculationPhiSum implements PotentialCalculationMolecular {
	 protected IVectorMutable fieldE;
	 protected final IVectorMutable ei,ej;
	 protected IVectorMutable Ai;
	 protected IVectorMutable Aj;
	 protected double secondDerivativeSum= 0;
	 protected DipoleSource dipoleSource;

	public PotentialCalculationPhiSum(ISpace space) {
		fieldE = space.makeVector();
		Ai = space.makeVector();
	    Aj = space.makeVector();
	    ei = space.makeVector();
	    ej = space.makeVector();
	}

	@Override
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		// TODO Auto-generated method stub 
		
		
		
		
		
		
		
		//TODO very important should come back in the future.
		
	}

	@Override
	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if(!(potential instanceof IPotentialMolecularSecondDerivative)){
			return;
		}
		
		IPotentialMolecularSecondDerivative potentialSeconDerivative = (IPotentialMolecularSecondDerivative) potential;
		
		Tensor[] t = potentialSeconDerivative.secondDerivative(molecules);
		
		IMolecule molecule0 = molecules.getMolecule(0);
		IMolecule molecule1 = molecules.getMolecule(1);
		ei.E(dipoleSource.getDipole(molecule0));
		ej.E(dipoleSource.getDipole(molecule1));
		ei.normalize();
		ej.normalize();
		
		
//		debug only  
//		IVectorMutable pos0 = atom1.getPosition();
//		IVectorMutable pos1 = atom0.getPosition();
//		Ai.E(pos1);
//		Ai.ME(pos0);
//		System.out.println("r = " + Ai);
//		System.out.println("ei = " + ei);
//		System.out.println("ej = " + ej);
		
		
		double si,sj;
		for(int k = 0; k < 3;k++){
//		for(int k = 2; k < 3;k++){ //debug only 						TODO!!!!!!!!!!!!!!>>>!>!>!>!>!!>!
			fieldE.E(0);
			fieldE.setX(k, 1);
			
			//phij*(1-xi^2)(1-xj^2) =  d2udthetaidthetaj.(fieldE cross ej).(field cross ei)*sj*si/sj/si
			
			
			Ai.E(fieldE);
			Ai.XE(ei);
			Aj.E(fieldE);
			Aj.XE(ej);
//			si = Math.sqrt(1-ei.getX(k)*ei.getX(k));
//			sj = Math.sqrt(1-ej.getX(k)*ej.getX(k));
//			System.out.println("fieldE =" + fieldE);
//			System.out.println("Ai =" + Ai);
//			System.out.println("Aj =" + Aj);
			
			t[0].transform(Aj);
			secondDerivativeSum += 2.0*Aj.dot(Ai);		//ij 
			
			
//			System.out.println("k = "+ k + " Phij = " + Ai.dot(Aj));//debug only for phij
			
//			debug only
//			System.out.println("dudij = \n" + t[0]);
//			System.out.println("dudij*Aj = " + Aj);
//			System.out.println("dudij*Aj*Ai = " + Aj.dot(Ai));

		//	debug only
//			Aj.E(fieldE);
//			Aj.XE(ej);
//			t[0].transpose();
//			t[0].transform(Ai);
//			secondDerivativeSum += Ai.dot(Aj);//ji
//			System.out.println("dudji = \n" + t[0]);
//			System.out.println("dudji*Ai = " + Ai);
//			System.out.println("dudji*Ai*Aj = " + Aj.dot(Ai));
//			System.exit(2);
			
			Ai.E(fieldE);
			Ai.XE(ei);
			Aj.E(Ai);
			t[1].transform(Ai);
			secondDerivativeSum += Ai.dot(Aj);		//ii
			
//			System.out.println("k = "+ k + " Phii = " + Ai.dot(Aj)/(si*si*si*si));  //debug only
			
			
			
			//debug only
//			System.out.println("dudii = \n" + t[1]);
//			System.out.println("dudii*Ai = " + Ai);
//			System.out.println("dudii*Ai*Ai = " + Ai.dot(Aj));
			
			
			Aj.E(fieldE);
			Aj.XE(ej);
			Ai.E(Aj);
			t[2].transform(Ai);
			secondDerivativeSum += Ai.dot(Aj);		//jj
			
//			System.out.println("k = "+ k + " Phjj = " + Ai.dot(Aj)/(sj*sj*sj*sj)); //bebug only for phjj
			
			//debug only
//			System.out.println("dudjj = \n" + t[2]);
//			System.out.println("dudjj*Aj = " + Ai);
//			System.out.println("dudjj*Aj*Aj = " + Ai.dot(Aj));
		}
		
//		System.exit(2);//TODO debug only
		
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

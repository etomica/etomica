package etomica.normalmode;

import etomica.api.*;
import etomica.box.Box;
import etomica.meam.PotentialCuLREP;
import etomica.meam.PotentialEAM;
import etomica.meam.PotentialEAM_LS;
import etomica.meam.PotentialEFS;
import etomica.potential.PotentialCalculation;
import etomica.space.ISpace;

public class PotentialCalculationEFSSP implements PotentialCalculation {
		
	public PotentialCalculationEFSSP(ISpace space, Box box, CoordinateDefinition coordinateDefinition, double temperature, double f1, boolean isLS){
		sum = new double[1];
//		sum = new double[7];
        this.temperature = temperature;
        rij = space.makeVector();
        Rij = space.makeVector();
        dri = space.makeVector();
        drj = space.makeVector();
        drij = space.makeVector();
        T1ij = space.makeVector();
        T1ik = space.makeVector();
        this.f1 = f1;
        this.isLS = isLS;
        this.box = box;
        this.coordinateDefinition = coordinateDefinition;
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().getAtomCount();
	}

	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        int nNbrAtoms = atoms.getAtomCount();
		IVector[] g = null;
		if(isLS){
			PotentialEAM_LS potentialSoft = (PotentialEAM_LS)potential;
	        g = potentialSoft.gradient(atoms);// gradient do nearestImage() !!!
		}else{
	        if(potential instanceof PotentialEAM){
	        	PotentialEAM potentialSoft = (PotentialEAM)potential;
		        g = potentialSoft.gradient(atoms);// gradient do nearestImage() !!!
	        }else if(potential instanceof PotentialEFS){
				PotentialEFS potentialSoft = (PotentialEFS)potential;
		        g = potentialSoft.gradient(atoms);// gradient do nearestImage() !!!
	       
	        }else if(potential instanceof PotentialCuLREP){
				PotentialCuLREP potentialSoft = (PotentialCuLREP)potential;
		        g = potentialSoft.gradient(atoms);// gradient do nearestImage() !!!
	        }
			
		}
		IVector ri = atoms.getAtom(0).getPosition();
		IVector Ri = coordinateDefinition.getLatticePosition(atoms.getAtom(0));
        dri.Ev1Mv2(ri, Ri);

        for (int j=0;j<nNbrAtoms;j++){//START from "1" NOT "0" because we need j != i
        	IVector rj = atoms.getAtom(j).getPosition();
        	IVector Rj = coordinateDefinition.getLatticePosition(atoms.getAtom(j));
        	rij.Ev1Mv2(ri , rj);
        	box.getBoundary().nearestImage(rij);
        	Rij.Ev1Mv2(Ri , Rj);
        	box.getBoundary().nearestImage(Rij);
        	drj.Ev1Mv2(rj , Rj);

        	drij.Ev1Mv2(dri , drj);
	        T1ij.Ea1Tv1(1.0/3.0/volume, Rij);
	        T1ij.PEa1Tv1(f1, drij); 
           	if(j != 0){
    	        sum[0] += g[j].dot(drij); //fij.drij
        	}//j!=0 
        }//j atoms loop
	}

	protected double function(IVector rij, IVector T1ij, IVector T2ij, double dW, double d2W){
		double rij2 = rij.squared();
		double rDr = dW/rij2 * T1ij.dot(T2ij) + (d2W-dW)/rij2/rij2 * T1ij.dot(rij)*(T2ij.dot(rij)); //B_new!!        		
		return rDr;
	}
	
	/**
	 * Sets the virial sum to zero, typically to begin a new virial-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public void reset() {
		for(int i=0;i<sum.length;i++){
			sum[i] = 0.0;
		}
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[] getSum() {
//		System.out.println(sum[2] +"  ,   "+ sum[3]+">>>>>"+(sum[2]+sum[3]));
		return sum;
	}
	
	private double[] sum;
	protected double volume , temperature;
	protected double f1;
	protected int nMol;
    protected final IVectorMutable rij;
    protected final IVectorMutable Rij;
    protected final IVectorMutable drj;
    protected final IVectorMutable dri;
    protected final IVectorMutable drij;
    protected final IVectorMutable T1ij, T1ik;

    protected final boolean isLS;

    protected final Box box;
    protected final CoordinateDefinition coordinateDefinition;
 }//end VirialSum

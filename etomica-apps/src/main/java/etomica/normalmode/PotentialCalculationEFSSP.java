/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.meam.PotentialCuLREP;
import etomica.meam.PotentialEAM;
import etomica.meam.PotentialEAM_LS;
import etomica.meam.PotentialEFS;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.space.Vector;

public class PotentialCalculationEFSSP implements PotentialCalculation {
		
	public PotentialCalculationEFSSP(Space space, Box box, CoordinateDefinition coordinateDefinition, double temperature, double f1, boolean isLS){
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
        nMol = coordinateDefinition.getBox().getLeafList().size();
	}

	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        int nNbrAtoms = atoms.size();
		Vector[] g = null;
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
		Vector ri = atoms.get(0).getPosition();
		Vector Ri = coordinateDefinition.getLatticePosition(atoms.get(0));
        dri.Ev1Mv2(ri, Ri);

        for (int j=0;j<nNbrAtoms;j++){//START from "1" NOT "0" because we need j != i
        	Vector rj = atoms.get(j).getPosition();
        	Vector Rj = coordinateDefinition.getLatticePosition(atoms.get(j));
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

	protected double function(Vector rij, Vector T1ij, Vector T2ij, double dW, double d2W){
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
    protected final Vector rij;
    protected final Vector Rij;
    protected final Vector drj;
    protected final Vector dri;
    protected final Vector drij;
    protected final Vector T1ij, T1ik;

    protected final boolean isLS;

    protected final Box box;
    protected final CoordinateDefinition coordinateDefinition;
 }//end VirialSum

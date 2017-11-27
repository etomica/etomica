/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.space.Vector;

public class PotentialCalculationLJSP implements PotentialCalculation {
		
	public PotentialCalculationLJSP(Space space, Box box, CoordinateDefinition coordinateDefinition, double temperature, double dP, double f1, double fe, double fee){
		sum = new double[1];
//		sum = new double[14];
        this.space = space;
        this.temperature = temperature;
        rij = space.makeVector();
        Rij = space.makeVector();
        drj = space.makeVector();
        dri = space.makeVector();
        drij = space.makeVector();
        drj_tmp = space.makeVector();
        T1ij = space.makeVector();
        T2ij = space.makeVector();
        T3ij = space.makeVector();
        this.dP = dP;
        this.f1 = f1;
        this.fe  = fe;
        this.fee = fee;
        
        this.box = box;
        this.coordinateDefinition = coordinateDefinition;
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().getAtomCount();
	}
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		Potential2SoftSpherical potentialSoft = (Potential2SoftSpherical)potential;
        Vector[] g = potentialSoft.gradient(atoms);// gradient do nearestImage() !!!
        int nNbrAtoms = atoms.getAtomCount();
        Vector ri = atoms.getAtom(0).getPosition();
        Vector Ri = coordinateDefinition.getLatticePosition(atoms.getAtom(0));
        dri.Ev1Mv2(ri, Ri);
        Vector rj;

        
        for (int j=1;j<nNbrAtoms;j++){//START from "1" NOT "0" because we need j != i
        	rj = atoms.getAtom(j).getPosition();
        	Vector Rj = coordinateDefinition.getLatticePosition(atoms.getAtom(j));
        	Rij.Ev1Mv2(Ri , Rj);
        	Vector shift_Rij = box.getBoundary().centralImage(Rij);
        	box.getBoundary().nearestImage(Rij);
        	rij.Ev1Mv2(ri , rj);
        	rij.PE(shift_Rij);
        	
//			double rij2 = rij.squared();
        	
        	drj.Ev1Mv2(rj , Rj);
        	drij.Ev1Mv2(dri , drj);
        	
//	        T1ij.Ea1Tv1(1.0/3.0/volume, Rij);
//	        T1ij.PEa1Tv1(f1, drij); 
//	        
//	        T2ij.setX(0, rij.getX(0));; 
//	        T2ij.setX(1, 0); 
//	        T2ij.setX(2, 0); 
//
//	        
//	        T3ij.setX(0, fe*drij.getX(0) + rij.getX(0));; 
//	        T3ij.setX(1, fe*drij.getX(1)); 
//	        T3ij.setX(2, fe*drij.getX(2)); 


//			double dW  = potentialSoft.du(rij2);
//	        double d2W = potentialSoft.d2u(rij2);
	        sum[0] += g[j].dot(drij); //fij.drij
//	        sum[1] += g[j].dot(rij); //fij.rij
//	        sum[2] += g[j].dot(Rij); //fij.Rij
//	        
//	        sum[3] += function(rij, T1ij, T1ij, dW, d2W); //Br
//	        sum[4] += function(rij, T1ij, drij, dW, d2W); //(dP/dT)r
//	        sum[5] += function(rij, rij,  rij,  dW, d2W); //B_Fluc
//	        
////Cij	
//	        sum[11] += function(rij, T2ij, T2ij, dW, d2W); //Br
//	        sum[12] += function(rij, T3ij, T3ij, dW, d2W); //Br
//	        sum[13] += g[j].getX(0) * drij.getX(0);//fxdrx:
//
//			DLJ.Ev1v2(rij,rij);
//	        DLJ.TE(1.0/(rij2*rij2)*(dW -  d2W));
//	        DLJ.PEa1Tt1(-dW/rij2,identity);
//	        drj_tmp.E(drj);
//	        DLJ.transform(drj_tmp); 
//	        double driDdrj = dri.dot(drj_tmp);
//	        double drjDdrj = drj.dot(drj_tmp);
//	        drj_tmp.E(dri);
//	        DLJ.transform(drj_tmp);
//	        double driDdri = dri.dot(drj_tmp);
//	        sum[6] += (driDdrj - 0.5*driDdri  - 0.5*drjDdrj); // Cv_r
//	        
////Sigmaij
//	        //xx
//	        sum[7] += g[j].getX(0)* rij.getX(0);//fxrx:
//	        sum[8] += g[j].getX(0)* Rij.getX(0);//fxRx:
//	        //xz
//	        sum[9]  +=   g[j].getX(0)*rij.getX(2);//fxrz
//	        sum[10] +=   g[j].getX(0)*Rij.getX(2);//fxRz
//	        
	        

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
	protected double dP,  f1, fe, fee;
	protected int nMol;
    protected final Vector rij;
    protected final Vector Rij;
    protected final Vector drj;
    protected final Vector dri ;
    protected final Vector drij;
    protected final Vector drj_tmp;
    protected final Vector T1ij, T2ij, T3ij;

    

    protected final Space space;
    protected final Box box;
    protected final CoordinateDefinition coordinateDefinition;
 }//end VirialSum

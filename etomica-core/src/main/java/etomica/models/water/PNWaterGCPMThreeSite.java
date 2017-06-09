/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * GCPM Water potential with three association site for Wertheim approach.  
 * 
 * @author Hye Min Kim
 */
public class PNWaterGCPMThreeSite extends PNWaterGCPM {
	
	private double mincosThetaH;
	private int bondType;
	public static final int BONDAC = 1, BONDBC = 2, BONDCA = 3, BONDCB = 4;
	protected boolean isAssociation;
	protected boolean bondingAngleRestriction;
	double minE;

	public PNWaterGCPMThreeSite(Space space, double minE, boolean isAssociation, boolean bondingAngleRestriction) {
		super(space);
		this.isAssociation = isAssociation;
		this.bondingAngleRestriction = bondingAngleRestriction;
		this.minE = minE;
	}
	
	public PNWaterGCPMThreeSite(Space space, double minE, boolean isAssociation) {
		super(space);
		this.isAssociation = isAssociation;
		this.minE = minE;
	}
	
	public double energy(IMoleculeList atoms){
        IAtomList water1Atoms = atoms.getMolecule(0).getChildList();
        IAtomList water2Atoms = atoms.getMolecule(1).getChildList();

        Vector O1r = water1Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        Vector O2r = water2Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        
        work.Ev1Mv2(O1r, O2r);
        
        double r2 = work.squared();
        
        if(r2<=core || r2 >= 12.25) {//outer shell = 3.5A
        	if(isAssociation) return 0.0;
        	else return super.energy(atoms);
        } 
        
        work.normalize();

        Vector H11r = water1Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        Vector H12r = water1Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();
        Vector H21r = water2Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        Vector H22r = water2Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();

        Vector M1r = water1Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        Vector M2r = water2Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        
        Vector rH11O2 = space.makeVector();
		rH11O2.E(H11r);
		rH11O2.ME(O2r);
		double distanceH11O2 = Math.sqrt(rH11O2.squared());
		
		Vector rH12O2 = space.makeVector();
		rH12O2.E(H12r);
		rH12O2.ME(O2r);
		double distanceH12O2 = Math.sqrt(rH12O2.squared());
		
		Vector rO1H21 = space.makeVector();
		rO1H21.E(O1r);
		rO1H21.ME(H21r);
		double distanceO1H21 = Math.sqrt(rO1H21.squared());
		
		Vector rO1H22 = space.makeVector();
		rO1H22.E(O1r);
		rO1H22.ME(H22r);
		double distanceO1H22 = Math.sqrt(rO1H22.squared());
		
        double[] distance = {distanceH11O2, distanceH12O2, distanceO1H21, distanceO1H22};
        double minDistance = distance[0];
        int minBond = 0;
        
		double cosTheta1=0;
		double cosTheta3=0;		
        
        for (int i = 1; i < distance.length; i++){
        	if (distance[i] < minDistance){
        		minDistance = distance[i];
        		minBond = i;
        	}
        }
        
		Vector O1H11 = space.makeVector();
		Vector H11O1 = space.makeVector();
		Vector O1H12 = space.makeVector();
		Vector H12O1 = space.makeVector();
		Vector H11O2 = space.makeVector();
		Vector O2H11 = space.makeVector();
		Vector H12O2 = space.makeVector();
		Vector O2H12 = space.makeVector();
		Vector O2H21 = space.makeVector();
		Vector H21O2 = space.makeVector();
		Vector O1H1A = space.makeVector();
		Vector O1H1 = space.makeVector();
		Vector O2H2A = space.makeVector();
		Vector O2H2 = space.makeVector();
		Vector H21O1 = space.makeVector();
		Vector O1H21 = space.makeVector();
		Vector H22O1 = space.makeVector();
		Vector O1H22 = space.makeVector();
		Vector O2H22 = space.makeVector();
		Vector H22O2 = space.makeVector();
		Vector O1O2 = space.makeVector();
		Vector O2O1 = space.makeVector();
		Vector O1OPrime = space.makeVector();
		Vector O2OPrime = space.makeVector();
		Vector O2H11Prime = space.makeVector();
		Vector O2H12Prime = space.makeVector();
		Vector O1H21Prime = space.makeVector();
		Vector O1H22Prime = space.makeVector();
		Vector crossProduct1 = space.makeVector();
		Vector crossProduct2 = space.makeVector();
		Vector projectionO1 = space.makeVector();
		Vector projectionO2 = space.makeVector();
		Vector projectionH11 = space.makeVector();
		Vector projectionH12 = space.makeVector();
		Vector projectionH21 = space.makeVector();
		Vector projectionH22 = space.makeVector();
		O1H11.E(H11r);
		O1H11.ME(O1r);
		O1H11.normalize();
		H11O1.E(O1r);
		H11O1.ME(H11r);
		H11O1.normalize();
		O1H12.E(H12r);
		O1H12.ME(O1r);
		O1H12.normalize();
		H12O1.E(O1r);
		H12O1.ME(H12r);
		H12O1.normalize();
		
		crossProduct1.E(O1H11);
		crossProduct1.XE(O1H12);
		crossProduct1.normalize();
		
		O1H1A.E(O1H11);
		O1H1A.PE(O1H12);
		O1H1A.TE(-1);
		O1H1A.normalize();
		
		O1H1.E(O1H11);
		O1H1.PE(O1H12);
		O1H1.normalize();
		
		O1H21.E(H21r);
		O1H21.ME(O1r);
		O1H21.normalize();
		H21O1.E(O1r);
		H21O1.ME(H21r);
		H21O1.normalize();
		O1H22.E(H22r);
		O1H22.ME(O1r);
		O1H22.normalize();
		H22O1.E(O1r);
		H22O1.ME(H22r);
		H22O1.normalize();
		O2H11.E(H11r);
		O2H11.ME(O2r);
		O2H11.normalize();
		H11O2.E(O2r);
		H11O2.ME(H11r);
		H11O2.normalize();
		O2H12.E(H12r);
		O2H12.ME(O2r);
		O2H12.normalize();
		H12O2.E(O2r);
		H12O2.ME(H12r);
		H12O2.normalize();
		O1O2.E(O2r);
		O1O2.ME(O1r);
		O1O2.normalize();
		O2O1.Ea1Tv1(-1, O1O2);
		O2H21.E(H21r);
		O2H21.ME(O2r);
		O2H21.normalize();
		H21O2.E(O2r);
		H21O2.ME(H21r);
		H21O2.normalize();
		O2H22.E(H22r);
		O2H22.ME(O2r);
		O2H22.normalize();
		H22O2.E(O2r);
		H22O2.ME(H22r);
		H22O2.normalize();
		
		crossProduct2.E(O2H21);
		crossProduct2.XE(O2H22);
		crossProduct2.normalize();
		
		O2H2A.E(O2H21);
		O2H2A.PE(O2H22);
		O2H2A.TE(-1);
		O2H2A.normalize();
		
		O2H2.E(O2H21);
		O2H2.PE(O2H22);
		O2H2.normalize();
		
		double fractionO1O2 = crossProduct1.dot(O1O2);
		double fractionO2O1 = crossProduct2.dot(O2O1);
		double fractionO1H21 = crossProduct1.dot(O1H21);
		double fractionO1H22 = crossProduct1.dot(O1H22);
		double fractionO2H11 = crossProduct2.dot(O2H11);
		double fractionO2H12 = crossProduct2.dot(O2H12);
		projectionO1.Ea1Tv1(fractionO2O1,crossProduct2); 
		projectionO2.Ea1Tv1(fractionO1O2, crossProduct1);
		projectionH11.Ea1Tv1(fractionO2H11, crossProduct2);
		projectionH12.Ea1Tv1(fractionO2H12, crossProduct2);
		projectionH21.Ea1Tv1(fractionO1H21, crossProduct1);
		projectionH22.Ea1Tv1(fractionO1H22, crossProduct1);
		O1OPrime.Ev1Mv2(O1O2, projectionO2);
		O2OPrime.Ev1Mv2(O2O1, projectionO1);
		O1H21Prime.Ev1Mv2(O1H21,projectionH21);
		O1H22Prime.Ev1Mv2(O1H22,projectionH22);
		O2H11Prime.Ev1Mv2(O2H11,projectionH11);
		O2H12Prime.Ev1Mv2(O2H12,projectionH12);


		
		if (minBond == 0){
        	cosTheta1 = H11O1.dot(H11O2);
        	cosTheta3 = O2H11.dot(O2H2);
		}
		if (minBond == 1){
        	cosTheta1 = H12O1.dot(H12O2);
        	cosTheta3 = O2H12.dot(O2H2);
		}
		if (minBond == 2){
        	cosTheta1 = H21O2.dot(H21O1);
        	cosTheta3 = O1H21.dot(O1H1);
		}
		if (minBond == 3){
        	cosTheta1 = H22O2.dot(H22O1);;
        	cosTheta3 = O1H22.dot(O1H1);
		}
		
        if (isAssociation&&bondType !=(minBond+1)){
        	return 0.0;
        }
        if (bondingAngleRestriction && isAssociation && (cosTheta1>-0.5||cosTheta3>0)){//bonding angle criteria
        	return 0.0;
        }


        
        double energy = super.energy(atoms); 
        

        if (isAssociation){
        	if (energy > -minE){
        		return 0.0;
        	}
        }
        else if(energy < -minE){
			return 0.0;
		}
 
		return energy;
	}
        
        public double getTheta() {return Math.acos(mincosThetaH);}
        
        /**
         * Accessor method for angle (in radians) describing width of cone.
         */
        public void setTheta(double t) {
            mincosThetaH   = Math.cos(t);
        }
        
        public void setBondType(int a){
        	bondType = a;
        }
}

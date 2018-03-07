/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;

public class P2HardAssociationGCPMReference extends PNWaterGCPM {
	
	private double mincosThetaH;
	private int bondType;
	public static final int BONDAC = 1, BONDBC = 2, BONDCA = 3, BONDCB = 4;
	protected boolean isAssociation;

	public P2HardAssociationGCPMReference(Space space, boolean isAssociation) {
		super(space);
		this.isAssociation = isAssociation;
	}
	
	public double energy(IMoleculeList atoms){
        IAtomList water1Atoms = atoms.get(0).getChildList();
        IAtomList water2Atoms = atoms.get(1).getChildList();

        Vector O1r = water1Atoms.get(SpeciesWater4P.indexO).getPosition();
        Vector O2r = water2Atoms.get(SpeciesWater4P.indexO).getPosition();
        
        work.Ev1Mv2(O1r, O2r);
        
        double r2 = work.squared();
        double eTot = 0.0;
      
        work.normalize();

        Vector H11r = water1Atoms.get(SpeciesWater4P.indexH1).getPosition();
        Vector H12r = water1Atoms.get(SpeciesWater4P.indexH2).getPosition();
        Vector H21r = water2Atoms.get(SpeciesWater4P.indexH1).getPosition();
        Vector H22r = water2Atoms.get(SpeciesWater4P.indexH2).getPosition();

        Vector M1r = water1Atoms.get(SpeciesWater4P.indexM).getPosition();
        Vector M2r = water2Atoms.get(SpeciesWater4P.indexM).getPosition();
        
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
        
        for (int i = 1; i < distance.length; i++){
        	if (distance[i] < minDistance){
        		minDistance = distance[i];
        		minBond = i;
        	}
        }

        if (r2 > core && r2 <12.25 &&isAssociation && bondType == (minBond+1)){//outer shell = 3.5A
        	eTot = Double.POSITIVE_INFINITY;
        }
 
		return eTot;
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

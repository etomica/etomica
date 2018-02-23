/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Degree;

public class PRotConstraint extends PotentialMolecular{

	/**
	 *  This class is created to constraint the rotation of the molecule
	 *  Make a hard-wall at 90 deg from the molecule's nominal orientation
	 *  
	 *  
	 */
	private static final long serialVersionUID = 1L;
	
	public PRotConstraint(Space space, CoordinateDefinitionNitrogen coordinateDefinition, Box box) {
		super(1, space);
		this.box = box;
		int numMolec = box.getMoleculeList().getMoleculeCount();
		
		molecOrientation = space.makeVector();
		initMolecOrientation = new Vector[numMolec][3];
		/*
		 * initializing the initial orientation of the molecule
		 */
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDefinition.getMoleculeOrientation(box.getMoleculeList().getMolecule(i));
		}
	
		
	}//End of Constructor

	
	public double energy(IMoleculeList molecules) {
						
		if(!onSwitch){
			return 0.0;
		}
		
		IMolecule molecule = molecules.getMolecule(0);
		int index = molecule.getIndex();
				
		Vector leafPos0 = molecule.getChildList().get(0).getPosition();
		Vector leaftPos1 = molecule.getChildList().get(1).getPosition();
		
		molecOrientation.Ev1Mv2(leaftPos1, leafPos0);
		molecOrientation.normalize();
		
		double cosangle = molecOrientation.dot(initMolecOrientation[index][0]);
		
		// to take care of machine precision issue
		// if cosangle = 1.0000000008; just set it to 1.0
		if (cosangle > 1.0 && cosangle < 1.0+1e-8){ 
			cosangle = 1.0;
		} 
		
		double angle = Math.acos(cosangle);
		
		if(  Degree.UNIT.fromSim(angle)<= constraintAngle ){
			return 0.0;
		
		} else {
			++ counter;
			return Double.POSITIVE_INFINITY;
		}
	}


	public double getRange() {
		return 0;
	}


	public void setBox(Box box) {
		this.box = box;
		
	}
	
	public void setConstraintAngle(double angle){
		this.constraintAngle = angle;
	}

	public void setSwitch(boolean switchOn){
		onSwitch = switchOn; 
	}
	
	
	private Vector[][] initMolecOrientation;
	private Vector molecOrientation;
	private Box box;
	protected double constraintAngle = 90.0; //in degree
	protected int counter=0;
	protected boolean onSwitch = true;

	
}

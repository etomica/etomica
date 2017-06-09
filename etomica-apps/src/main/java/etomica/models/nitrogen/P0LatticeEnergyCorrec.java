/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMolecular;
import etomica.space.Space;
import etomica.units.Kelvin;

/**
 *  The lattice energy correction for Nitrogen model
 *  zero-body interaction
 *  
 *  All the parameters are determined from the fitted equation. 
 *   (uLattice_infinite - uLattice) as function of rho
 * 
 * @author Tai Boon Tan
 *
 */
public class P0LatticeEnergyCorrec extends PotentialMolecular{

	public P0LatticeEnergyCorrec(Space space){
		super(0, space);
		coeff = new double[3];
	}
	
	public double energy(IMoleculeList atoms) {
		double rho = numMolec/box.getBoundary().volume();
		return uCorrection(rho);
	}
	
	public void setBox(Box p){
		this.box = p;
	   if(species==null){
       	throw new RuntimeException("<P0LatticeEnergyCorrec.java> Must set Species First");
       }
                      
       numMolec = p.getNMolecules(species);
       
       if (numMolec == 32){
       	caseNumMolec = 1;
       	coeff[0] =  0.288064;
       	coeff[1] = -8.03294;
       	coeff[2] = - 490092;
       	
       } else if (numMolec == 108){
       	caseNumMolec = 2;
       	coeff[0] = -1.14609;
       	coeff[1] =  226.723;
       	coeff[2] = - 153683;
       	
       } else if (numMolec == 256){
       	caseNumMolec = 3;
       	coeff[0] = -0.210368;
       	coeff[1] =   34.4197;
       	coeff[2] = - 74137.8;
           
       }else if (numMolec == 500){
       	caseNumMolec = 4;
       	coeff[0] = -0.0761822;
       	coeff[1] =    15.1892;
       	coeff[2] = -  31121.6;
           
       } 
	}
	
    private double uCorrection(double rho){
    	double rho2 = rho * rho;
    	
    	if (caseNumMolec==0){
    		return  0.0;
    	
    	} else {
    		/*
    		 * return the correction energy for the total system
    		 * NOT the correction energy per molecule
    		 * The coeff was fitted with energy in K
    		 * so we have to convert the unit to simulation unit 
    		 */
    		return  Kelvin.UNIT.toSim(numMolec*(coeff[0] + coeff[1]*rho + coeff[2]*rho2));
    	
    	}
    }
	
	public ISpecies getSpecies() {
		return species;
	}

	public void setSpecies(ISpecies species) {
		this.species = species;
	}
	
	public double getRange() {
		return 0;
	}

	private ISpecies species;
	private Box box;
	private double[] coeff;
    private int caseNumMolec = 0;
    protected int numMolec;
	private static final long serialVersionUID = 1L;


}

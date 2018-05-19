/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.traPPE;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for TraPPE model of methanol: http://www.chem.umn.edu/groups/siepmann/trappe/molname.php#
 * 
 * CH3 group of methanol is united into a single site 
 * 
 * K.R. Schadel 2008
 */
public class ConformationMethanol implements IConformation {

    public ConformationMethanol(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
    	
    	double bondCH3O = 1.43; // Angstroms
    	double bondOH = 0.945; // Angstroms  (Chen et al report 0.945 Angstroms..., the website says 0.95 Angstroms)
    	double angleEq = 108.50*Math.PI/180; // equilibrium bond angle in radians (mcWiggle will change this appropriately)
    	
    	IAtom cH3 = list.get(SpeciesMethanol.indexCH3);
        cH3.getPosition().E(new double[] {bondCH3O, 0.0, 0.0});
        
        IAtom oxygen = list.get(SpeciesMethanol.indexO);
        oxygen.getPosition().E(new double[] {0.0, 0.0, 0.0});
        
    	// hydrogen attached to oxygen
        IAtom hydrogen = list.get(SpeciesMethanol.indexH);
        hydrogen.getPosition().E(new double[] {bondOH*Math.cos(angleEq), bondOH*Math.sin(angleEq), 0.0});
   
        
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;

}


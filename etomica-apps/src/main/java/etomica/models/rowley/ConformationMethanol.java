/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for methanol in standard orientation as published in Rowley et al (2006) 
 * 
 * OpenOffice Spreadsheet file, Methanol coordinates, used to calculate coordinates from bond lengths and angles: 
 * 
 * K.R. Schadel May 2008
 */
public class ConformationMethanol implements IConformation {

    public ConformationMethanol(Space space, boolean pointCharges) {
        this.space = space;
        this.pointCharges = pointCharges;
    }
    
    public void initializePositions(IAtomList list){
        
    	// hydrogen attached to oxygen
        IAtom alpha_hydrogen = list.get(SpeciesMethanol.indexaH);
        alpha_hydrogen.getPosition().E(new double[] {0.9114, 0.0, 1.7172});
        
        IAtom oxygen = list.get(SpeciesMethanol.indexO);
        oxygen.getPosition().E(new double[] {0.0, 0.0, 1.4202});
                
        IAtom carbon = list.get(SpeciesMethanol.indexaC);
        carbon.getPosition().E(new double[] {0.0, 0.0, 0.0});
        
        // hydrogen on opposite of molecule from alpha hydrogen
        IAtom hydrogen_1 = list.get(SpeciesMethanol.indexH1);
        hydrogen_1.getPosition().E(new double[] {-1.0393, 0.0, -0.3115});
        
        // hydrogen closer to alpha hydrogen
        IAtom hydrogen_2a = list.get(SpeciesMethanol.indexH2a);
        hydrogen_2a.getPosition().E(new double[] {0.4850, -0.8920, -0.4101});
        
        // hydrogen closer to alpha hydrogen
        IAtom hydrogen_2b = list.get(SpeciesMethanol.indexH2b);
        hydrogen_2b.getPosition().E(new double[] {0.4850,  0.8920, -0.4101});
        
        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        if(pointCharges) {
        	
        	// Methanol model with point charges at O, aC, aH
        	
        	// satellite site to model location of high electron density for oxygen
        	IAtom x = list.get(SpeciesMethanol.indexX);
        	x.getPosition().E(new double[] {-0.0753, 0.0, 1.4749});	
        } else {

        	// Methanol model without point charges
        	
        	// satellite site to model location of high electron density for oxygen
        	IAtom x = list.get(SpeciesMethanol.indexX);
        	x.getPosition().E(new double[] {-0.7774, 0.0, 1.9845});
        }
       
        
        
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected boolean pointCharges;

}


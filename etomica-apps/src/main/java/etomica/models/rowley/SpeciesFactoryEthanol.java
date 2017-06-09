/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.virial.SpeciesFactory;

/**
 * SpeciesFactory that makes ethanol.
 */
public class SpeciesFactoryEthanol implements SpeciesFactory, java.io.Serializable {
	
	public SpeciesFactoryEthanol(boolean pointCharges) {
		this.pointCharges = pointCharges;
		
	}
    
    public ISpecies makeSpecies(Space space) {
    	
    	// The satellite site, X, is closer to the oxygen atom in the model with point charges.
    	SpeciesEthanol species = new SpeciesEthanol(space, pointCharges);
        
        return species;
    }
    
    private static final long serialVersionUID = 1L;
    protected boolean pointCharges;

}


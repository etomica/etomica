/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.SpeciesSpheres;

/**
 * @author nsives, msellers
 * 
 * Sets carbon atom coordinates in tubular form.
 *  
 */

public class SpeciesTube extends SpeciesSpheres {
	
	
	public SpeciesTube(int atomsPerRing, int numberOfRings, Space _space){
		super(atomsPerRing*numberOfRings, new ElementSimple("T", Double.POSITIVE_INFINITY), 
                new ConformationTube(_space, atomsPerRing), _space);
		setIsDynamic(true);
	}
	
    private static final long serialVersionUID = 1L;
}

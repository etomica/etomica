/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on May 18, 2005
 */
package etomica.models.hexane;

import etomica.chem.elements.ElementSimple;
import etomica.space.Space;

/**
 * Species used to create hexane molecules per Dr. Monson's data. Hydrogen
 * molecules are ignored.
 * 
 * @author nancycribbin
 */

public class SpeciesHexane extends etomica.species.SpeciesSpheres {
    public SpeciesHexane(Space _space){
        super(6, new ElementSimple("M", 1.0), 
              new ConformationHexane(_space), _space);
        
        bondLength = 0.4;
//        bondAngle = 109.47 * 2.0 * Math.PI/360.0;
        phi = (180 - 109.47) / 360.0 * 2.0 * Math.PI;
        lowerLimit = 108.6919204 / 360.0 * 2.0 * Math.PI;
        upperLimit = 251.3080797 / 360.0 * 2.0 * Math.PI;
        
    }
    
    double bondLength, phi, upperLimit, lowerLimit;

    public double getBondAngle() {
        return phi;
    }

    public double getBondLength() {
        return bondLength;
    }

    public double getLowerLimit() {
        return lowerLimit;
    }

    public double getUpperLimit() {
        return upperLimit;
    }
}

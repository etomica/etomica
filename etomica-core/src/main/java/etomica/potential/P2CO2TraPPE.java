/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * I used this original EPM2model to describe TraPPECO2 model
 * the two models are different in parameters. For TraPPE, there is no k(theta) term.
 * the TraPPE CO2 is from Stubbs, Drake-Wilhelm and Siepmann 'Partial Molar Volume and Solvation structure of Na in SCCO2', 2005
 * 
 * @author shu
 * Oct.20,2010
 */
public class P2CO2TraPPE extends P2CO2EMP {

    public P2CO2TraPPE(Space space) {
    	
       //the following is EPM2 as comparison
    	//super(space, 2.757, 2.8921, 3.033, Kelvin.UNIT.toSim(28.129), Kelvin.UNIT.toSim(47.588), Kelvin.UNIT.toSim(80.507), Electron.UNIT.toSim(0.6512));
    	
        super(space, 2.8,2.925 , 3.05, Kelvin.UNIT.toSim(27), Kelvin.UNIT.toSim(46.1844), Kelvin.UNIT.toSim(79), Electron.UNIT.toSim(0.70));
       
    }

    private static final long serialVersionUID = 1L;
    
}

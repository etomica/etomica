/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.water;

import etomica.atom.MoleculePositionCOM;
import etomica.space.Space;
import etomica.units.Calorie;
import etomica.units.Electron;
import etomica.units.Mole;

/** 
 * TIP4P potential for water.  All the real work happens in P2Water4P.
 */
public class P2WaterTIP4P extends P2Water4P {

    // U = A/r^12 + C/r^6
    private static double A = 600e3; // kcal A^12 / mol
    private static double C = 610; // kcal A^6 / mol
    private static double s6 = A/C;
    public static double s = Math.pow(s6, 1.0/6.0);
    public static double e = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C/s6*1000))/4;
    public static double qH = Electron.UNIT.toSim(0.52);
    
    
    public P2WaterTIP4P(Space space) {
    	this(space, Double.POSITIVE_INFINITY);
    }
    public P2WaterTIP4P(Space space, double rCut) {
	    super(space, s, e, qH, rCut, new MoleculePositionCOM(space));
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

	public double getSigma() {return sigma;}		//TODO can I just add them here?

	public double getEpsilon() {return epsilon;}
	
	public double getChargeM() {return chargeM;}
	public double getChargeH() {return chargeH;}
	
	

    private static final long serialVersionUID = 1L;
}

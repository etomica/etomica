/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import etomica.modifier.Modifier;
import etomica.potential.P2SquareWell;
import etomica.units.Dimension;

/**
 * @author Matt Moynihan
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
class WellDepthModifier implements Modifier {

	P2SquareWell potential;

	WellDepthModifier(P2SquareWell pot) {
		potential = pot;
	}

	public String getLabel() {
		return "WellDepth";
	}

	public Dimension getDimension() {
		return etomica.units.Dimension.ENERGY;
	}

	public void setValue(double d) {
		potential.setEpsilon(d);
	}

	public double getValue() {
		return potential.getEpsilon();
	}
}

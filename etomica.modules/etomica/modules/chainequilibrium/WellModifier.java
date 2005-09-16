
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
class WellModifier implements Modifier {

	double currentValue;
	double fullDiameter;
	P2SquareWell potential;
	boolean initializing;
	
	WellModifier(P2SquareWell pot) {
		potential = pot;
		fullDiameter = potential.getCoreDiameter() * potential.getLambda();
//		 ********* Marker 
		System.out.println("Graphic: Created a Well Modifer");
	}

	public String getLabel() {
		return "WellMod";
	}

	public Dimension getDimension() {
		return etomica.units.Dimension.NULL;
	}

	public void setValue(double d) {
		if (initializing)
			return;
		currentValue = d;
		double x = 0.01 * currentValue;
		fullDiameter = potential.getCoreDiameter() * potential.getLambda();
		potential.setCoreDiameter(x * fullDiameter);
		potential.setLambda(1.0 / x);
	}

	public double getValue() {
//		 ********* Marker 
		System.out.println("Graphic: Well Modifer: get Value function");
		return currentValue;
	}
}//end of WellModulator

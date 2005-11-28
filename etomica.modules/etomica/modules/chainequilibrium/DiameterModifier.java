/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import etomica.atom.AtomTypeSphere;
import etomica.graphics.DisplayPhase;
import etomica.modifier.Modifier;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;

/**
 * @author Matt Moynihan
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
class DiameterModifier implements Modifier {
	P2SquareWellBonded potentialRR, potentialRB, potentialBB;
	SpeciesSpheresMono speciesR, speciesB;
	DisplayPhase display;

	DiameterModifier(P2SquareWellBonded potentialRR,
			P2SquareWellBonded potentialRB, P2SquareWellBonded potentialBB,
			SpeciesSpheresMono speciesR, SpeciesSpheresMono speciesB) {
		this.potentialRR = potentialRR;
		this.potentialRB = potentialRB;
		this.potentialBB = potentialBB;
		this.speciesR = speciesR;
		this.speciesB = speciesB;
	}

	public Dimension getDimension() {
		return etomica.units.Dimension.LENGTH;
	}

	public void setValue(double d) {
		if (d == 0.0)
			d = 0.01;
		double changeFraction = d
				/ (potentialRR.getCoreDiameter() * potentialRR.getLambda());
		double newCoreDiameter = changeFraction
				* potentialRR.getCoreDiameter();
		potentialRR.setCoreDiameter(newCoreDiameter);
		potentialRB.setCoreDiameter(newCoreDiameter);
		potentialBB.setCoreDiameter(newCoreDiameter);
		((AtomTypeSphere)speciesR.getMoleculeType()).setDiameter(d);
        ((AtomTypeSphere)speciesB.getMoleculeType()).setDiameter(d);
		if (display != null)
			display.repaint();
	}

	public double getValue() {
		return ((AtomTypeSphere)speciesR.getMoleculeType()).diameter(null);
	}

	public void setDisplay(DisplayPhase display) {
		this.display = display;
	}

	public String getLabel() {
		return "Diameter";
	}
}

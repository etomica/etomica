package etomica.chem.models;

import etomica.chem.electrostatics.Electrostatic;
import etomica.chem.elements.Element;
import etomica.potential.Potential;
import etomica.space.Space;

/**
 * @author zhaofang
 *
 */
public class Point extends ModelAtomic {

	/**
	 * 
	 */
	public Point() {
		super();
	}

	/**
	 * @param element
	 */
	public Point(Element element) {
		super(element);
	}

	/**
	 * @param element
	 * @param electrostatic
	 */
	public Point(Element element, Electrostatic electrostatic) {
		super(element, electrostatic);
	}

	
	/**
	 * @see etomica.chem.models.Model#makePotential(etomica.SimulationElement)
	 */
	public Potential makePotential(Space space) {
		return null;
	}

}

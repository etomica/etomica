package etomica.chem.models;

import etomica.chem.Electrostatic;
import etomica.chem.Element;
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
	 * @see etomica.chem.Model#makePotential(etomica.SimulationElement)
	 */
	public Potential makePotential(Space space) {
		return null;
	}

}

/*
 * Created on Jan 30, 2004
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package etomica.chem.models;

import etomica.Potential;
import etomica.SimulationElement;
import etomica.chem.Electrostatic;
import etomica.chem.Element;

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
	public Potential makePotential(SimulationElement parent) {
		return Potential.NullPotential((etomica.Simulation)parent);
	}

}

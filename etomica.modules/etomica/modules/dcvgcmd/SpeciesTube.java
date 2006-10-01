package etomica.modules.dcvgcmd;

import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomTypeLeaf;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.species.SpeciesSpheres;

/**
 * @author nsives, msellers
 * 
 * Sets carbon atom coordinates in tubular form.
 *  
 */

public class SpeciesTube extends SpeciesSpheres {
	
	
	public SpeciesTube(Simulation sim, int atomsPerRing, int numberOfRings){
		super(sim, atomsPerRing*numberOfRings, new ElementSimple("T", Double.POSITIVE_INFINITY), 
                new ConformationTube(sim.space, atomsPerRing));
	
		setNMolecules(1);
	}
	
	
}
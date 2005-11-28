package etomica.modules.dcvgcmd;

import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomTypeLeaf;
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
		super(sim);
	
		ConformationTube conformationTube = new ConformationTube(sim.space, atomsPerRing);
		
		factory.setConformation(conformationTube);
		
		((AtomTypeLeaf)factory.getType()).setMass(Double.POSITIVE_INFINITY);
	
		setNMolecules(1);
		((AtomFactoryHomo)factory).setAtomsPerGroup(atomsPerRing * numberOfRings);
	}
	
	
}
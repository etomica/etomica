package etomica.modules.dcvgcmd;

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
		
		setMass(Double.POSITIVE_INFINITY);
	
		setNMolecules(1);
	
		setAtomsPerMolecule(atomsPerRing * numberOfRings);
	}
	
	
}
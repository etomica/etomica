package etomica.modules.dcvgcmd;

import etomica.chem.elements.ElementSimple;
import etomica.simulation.ISimulation;
import etomica.species.SpeciesSpheres;

/**
 * @author nsives, msellers
 * 
 * Sets carbon atom coordinates in tubular form.
 *  
 */

public class SpeciesTube extends SpeciesSpheres {
	
	
	public SpeciesTube(ISimulation sim, int atomsPerRing, int numberOfRings){
		super(sim, atomsPerRing*numberOfRings, new ElementSimple("T", Double.POSITIVE_INFINITY), 
                new ConformationTube(sim.getSpace(), atomsPerRing));
	}
	
    private static final long serialVersionUID = 1L;
}

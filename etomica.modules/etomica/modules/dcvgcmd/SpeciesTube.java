package etomica.modules.dcvgcmd;

import etomica.api.ISimulation;

import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.SpeciesSpheres;

/**
 * @author nsives, msellers
 * 
 * Sets carbon atom coordinates in tubular form.
 *  
 */

public class SpeciesTube extends SpeciesSpheres {
	
	
	public SpeciesTube(ISimulation sim, int atomsPerRing, int numberOfRings, Space _space){
		super(sim, atomsPerRing*numberOfRings, new ElementSimple("T", Double.POSITIVE_INFINITY), 
                new ConformationTube(_space, atomsPerRing), _space);
	}
	
    private static final long serialVersionUID = 1L;
}

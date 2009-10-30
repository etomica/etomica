package etomica.modules.dcvgcmd;

import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheres;

/**
 * @author nsives, msellers
 * 
 * Sets carbon atom coordinates in tubular form.
 *  
 */

public class SpeciesTube extends SpeciesSpheres {
	
	
	public SpeciesTube(int atomsPerRing, int numberOfRings, ISpace _space){
		super(atomsPerRing*numberOfRings, new ElementSimple("T", Double.POSITIVE_INFINITY), 
                new ConformationTube(_space, atomsPerRing), _space);
		setIsDynamic(true);
	}
	
    private static final long serialVersionUID = 1L;
}

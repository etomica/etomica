/*
 * Created on Mar 29, 2005 
 */
package etomica.models.hexane;

import etomica.ConformationChainZigZag;
import etomica.Space;


/**
 * Defines the vectors used to create a hexane molecule according to Dr. Monson's
 * data
 * 
 * @author nancycribbin
 */
public class ConformationHexane extends ConformationChainZigZag {

	public ConformationHexane(Space space){
		super(space);
		
		 v1.setX(0, 0.188545552);
		 v1.setX(1, 0.296942436);
		 v1.setX(2, 0.1904619760);
		 
		 v2.setX(0, 0.037596211);
		 v2.setX(1, -0.092687873);
		 v2.setX(2, 0.387292503);
	}
}

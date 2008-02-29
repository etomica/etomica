/*
 * Created on Mar 24, 2005
 */
package etomica.config;

import etomica.api.IVector;
import etomica.space.Space;

/**
 * @author nancycribbin
 *
 */
public class ConformationChain2D extends ConformationChain {
	
	public ConformationChain2D(Space space, IVector[] vex){
		super(space);
		if(vex.length != vectors.length){
			throw new IllegalArgumentException("Different vector array lengths in ConformationChain2D.");
		}
				
		for(int i = 0; i < vex.length; i++){
			vectors[i].E(vex[i]);
		}
		tracker = 0;
	}
	
	/* (non-Javadoc)
	 * @see etomica.ConformationChain#reset()
	 */
	public void reset() { 
		tracker = 0;
	}

	/* (non-Javadoc)
	 * @see etomica.ConformationChain#nextVector()
	 */
	public IVector nextVector() {
		if(tracker<vectors.length){
			tracker += 1;
			return vectors[tracker-1];
		}
	    reset();
	    tracker += 1;
	    return vectors[tracker-1];
	}

	IVector[] vectors;
	int tracker;			//Tracker is used to track which vector the counter is on.

}

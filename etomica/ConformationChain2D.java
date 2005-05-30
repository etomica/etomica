/*
 * Created on Mar 24, 2005
 */
package etomica;

import etomica.space.Vector;

/**
 * @author nancycribbin
 *
 */
public class ConformationChain2D extends ConformationChain {
	
	public ConformationChain2D(Space space, Vector[] vex){
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
	public Vector nextVector() {
		if(tracker<vectors.length){
			tracker += 1;
			return vectors[tracker-1];
		}
	    reset();
	    tracker += 1;
	    return vectors[tracker-1];
	}

	Vector[] vectors;
	int tracker;			//Tracker is used to track which vector the counter is on.

}

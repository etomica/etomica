/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Sep 23, 2004 by kofke
 */
package etomica.potential;

import etomica.nbr.NeighborCriterion;

import java.util.Arrays;

/**
 * Instances of this class are created by AtomType classes to
 * store the list of potentials and criteria that apply to 
 * the atoms of that type.
 */
public class PotentialArray implements java.io.Serializable {

	public PotentialArray() {
		super();
	}

    /**
     * Returns the index that this type has assigned to the given
     * potential.  If potential has not been previously identified
     * to this instance, a new index will be returned.
     * @param potential The potential of interest
     * @return index associated with potential
     */
    public synchronized int getPotentialIndex(IPotential potential) {
        if(mostRecentIndex > -1 && potentials[mostRecentIndex]==potential) return mostRecentIndex;

        mostRecentIndex = -1;
        while(++mostRecentIndex < potentials.length) {
            if(potentials[mostRecentIndex]==potential) {
                return mostRecentIndex; 
            }
        }
        throw new IllegalArgumentException("Potential "+potential+" is unknown");
    }
    
    /**
     * Identifies the given potential to this type instance and returns
     * the index assigned to it.  If potential was previously added,
     * no action is taken other than to return the previously assigned
     * index for the potential (acts same as getPotentialIndex in this 
     * case).
     * @param newPotential the potential being added
     * @return the new or previously assigned index for the potential
     */
    public int addPotential(IPotential newPotential) {
    	for(mostRecentIndex=0; mostRecentIndex<potentials.length; mostRecentIndex++) {
    		if(potentials[mostRecentIndex] == newPotential) return mostRecentIndex;
    	}
        potentials = (IPotential[])etomica.util.Arrays.addObject(potentials, newPotential);
        // make room for a criterion to be added via setCriterion
        criteria = Arrays.copyOf(criteria, potentials.length);
    	return potentials.length-1;
    }
    
    /**
     * Sets the criterion associated with the given potential
     */
    public void setCriterion(IPotential potential, NeighborCriterion criterion) {
        criteria[getPotentialIndex(potential)] = criterion;
    }

    /**
     * Removes the given potential and the associated criterion from the list 
     * held by this type instance, and returns the index that was assigned to 
     * it.  If the potential was not previously added to this instance, returns
     * -1.
     */
    public int removePotential(IPotential potential) {
        for (int i=0; i<potentials.length; i++) {
    		if (potentials[i] == potential) {
    	    	Potential[] newPotentials = new Potential[potentials.length-1];
    	    	System.arraycopy(potentials,0,newPotentials,0,i);
    	    	System.arraycopy(potentials,i+1,newPotentials,i,potentials.length-i-1);
    	    	potentials = newPotentials;

                NeighborCriterion[] newCriteria = new NeighborCriterion[criteria.length-1];
                System.arraycopy(criteria,0,newCriteria,0,i);
                System.arraycopy(criteria,i+1,newCriteria,i,criteria.length-i-1);
                criteria = newCriteria;
    	    	return i;
    		}
    	}
    	return -1;
    }
    
    public final IPotential[] getPotentials() {
    	return potentials;
    }
    
    public final NeighborCriterion[] getCriteria() {
        return criteria;
    }
    
    private static final long serialVersionUID = 1L;
    private IPotential[] potentials = new IPotential[0];
    private NeighborCriterion[] criteria = new NeighborCriterion[0];
    private int mostRecentIndex = -1;

}

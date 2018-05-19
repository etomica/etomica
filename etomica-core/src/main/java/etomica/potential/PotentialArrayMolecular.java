/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Sep 23, 2004 by kofke
 */
package etomica.potential;

import etomica.nbr.molecule.NeighborCriterionMolecular;

import java.util.Arrays;

/**
 * Instances of this class are created by MoleculeType classes to
 * store the list of potentials and criteria that apply to 
 * the molecules of that type.
 */
public class PotentialArrayMolecular implements java.io.Serializable {

	public PotentialArrayMolecular() {
		super();
	}

    /**
     * Returns the index that this type has assigned to the given
     * potential.  If potential has not been previously identified
     * to this instance, a new index will be returned.
     * @param potential The potential of interest
     * @return index associated with potential
     */
    public int getPotentialIndex(IPotentialMolecular potential) {
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
    public int addPotential(IPotentialMolecular newPotential) {
    	for(mostRecentIndex=0; mostRecentIndex<potentials.length; mostRecentIndex++) {
    		if(potentials[mostRecentIndex] == newPotential) return mostRecentIndex;
    	}
        potentials = (IPotentialMolecular[])etomica.util.Arrays.addObject(potentials, newPotential);
        // make room for a criterion to be added via setCriterion
        criteria = Arrays.copyOf(criteria, potentials.length);
    	return potentials.length-1;
    }
    
    /**
     * Sets the criterion associated with the given potential
     */
    public void setCriterion(IPotentialMolecular potential, NeighborCriterionMolecular criterion) {
        criteria[getPotentialIndex(potential)] = criterion;
    }
    
    /**
     * Removes the given potential and the associated criterion from the list 
     * held by this type instance, and returns the index that was assigned to 
     * it.  If the potential was not previously added to this instance, returns
     * -1.
     */
    public int removePotential(IPotentialMolecular potential) {
        for (int i=0; i<potentials.length; i++) {
    		if (potentials[i] == potential) {
    	    	PotentialMolecular[] newPotentials = new PotentialMolecular[potentials.length-1];
    	    	System.arraycopy(potentials,0,newPotentials,0,i);
    	    	System.arraycopy(potentials,i+1,newPotentials,i,potentials.length-i-1);
    	    	potentials = newPotentials;

                NeighborCriterionMolecular[] newCriteria = new NeighborCriterionMolecular[criteria.length-1];
                System.arraycopy(criteria,0,newCriteria,0,i);
                System.arraycopy(criteria,i+1,newCriteria,i,criteria.length-i-1);
                criteria = newCriteria;
    	    	return i;
    		}
    	}
    	return -1;
    }
    
    public final IPotentialMolecular[] getPotentials() {
    	return potentials;
    }
    
    public final NeighborCriterionMolecular[] getCriteria() {
        return criteria;
    }
    
    private static final long serialVersionUID = 1L;
    private IPotentialMolecular[] potentials = new IPotentialMolecular[0];
    private NeighborCriterionMolecular[] criteria = new NeighborCriterionMolecular[0];
    private int mostRecentIndex = -1;

}

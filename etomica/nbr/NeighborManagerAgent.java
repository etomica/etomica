/*
 * History
 * Created on Sep 23, 2004 by kofke
 */
package etomica.nbr;

import etomica.Potential;

/**
 * Instances of this class are created by AtomType classes to
 * store the list of potentials and criteria that apply to 
 * the atoms of that type.
 */
public class NeighborManagerAgent {

	/**
	 * 
	 */
	public NeighborManagerAgent() {
		super();
	}

    /**
     * Returns the index that this type has assigned to the given
     * potential.  If potential has not been previously identified
     * to this instance, an ArrayIndexOutOfBoundsException is thrown.
     * This should be caught and then addPotential can be called
     * to identify the potential to this instance and obtain an index.
     * @param potential The potential of interest
     * @return index associated with potential
     */
    public int getPotentialIndex(Potential potential) {
    	if(potentials[mostRecentIndex]==potential) return mostRecentIndex;
    	else {
    		mostRecentIndex = -1;
    		while(true) {//potential not in list will lead to ArrayIndexOutOfBoundsException 
    			if(potentials[++mostRecentIndex]==potential) {
    				return mostRecentIndex; 
    			}
    		}
    	}
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
    public int addPotential(Potential newPotential) {
    	for(mostRecentIndex=0; mostRecentIndex<potentials.length; mostRecentIndex++) {
    		if(potentials[mostRecentIndex] == newPotential) return mostRecentIndex;
    	}
    	Potential[] newArray = new Potential[potentials.length+1];
    	System.arraycopy(potentials, 0, newArray, 0, potentials.length);
    	newArray[potentials.length] = newPotential;
    	potentials = newArray;
    	return potentials.length-1;
    }

    public int removePotential(Potential potential) {
    	for (int i=0; i<potentials.length; i++) {
    		if (potentials[i] == potential) {
    	    	Potential[] newPotentials = new Potential[potentials.length-1];
    	    	System.arraycopy(potentials,0,newPotentials,0,i);
    	    	System.arraycopy(potentials,i+1,newPotentials,i,potentials.length-i-1);
    	    	potentials = newPotentials;
    	    	return i;
    		}
    	}
    	return -1;
    }
    
    public Potential[] getPotentials() {
    	return potentials;
    }
    
    public NeighborCriterion getCriterion() {
    	return (criteria.length<1 ? null : criteria[0]);
    }
    public void addCriterion(NeighborCriterion criterion) {
    	for(int i=0; i<criteria.length; i++) {
    		if(criteria[i] == criterion) return;
    	}
    	NeighborCriterion[] newCriteria = new NeighborCriterion[criteria.length+1];
    	System.arraycopy(criteria,0,newCriteria,0,criteria.length);
    	newCriteria[criteria.length] = criterion;
    	criteria = newCriteria;
    }
    
    public boolean removeCriterion(NeighborCriterion criterion) {
    	for (int i=0; i<criteria.length; i++) {
    		if (criteria[i] == criterion) {
    	    	NeighborCriterion[] newCriteria = new NeighborCriterion[criteria.length-1];
    	    	System.arraycopy(criteria,0,newCriteria,0,i);
    	    	System.arraycopy(criteria,i+1,newCriteria,i,criteria.length-i-1);
    	    	criteria = newCriteria;
    	    	return true;
    		}
    	}
    	return false;
    }
    
    private Potential[] potentials = new Potential[0];
    private int mostRecentIndex = 0;
    private NeighborCriterion[] criteria = new NeighborCriterion[0];

}

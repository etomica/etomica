/*
 * History
 * Created on Sep 23, 2004 by kofke
 */
package etomica.potential;

import etomica.utility.Arrays;

/**
 * Instances of this class are created by AtomType classes to
 * store the list of potentials and criteria that apply to 
 * the atoms of that type.
 */
public class PotentialArray implements java.io.Serializable {

	/**
	 * 
	 */
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
    public int getPotentialIndex(Potential potential) {
        if(mostRecentIndex > -1 && potentials[mostRecentIndex]==potential) return mostRecentIndex;

        mostRecentIndex = -1;
        while(++mostRecentIndex < potentials.length) {
            if(potentials[mostRecentIndex]==potential) {
                return mostRecentIndex; 
            }
        }
        mostRecentIndex = addPotential(potential); 
        return mostRecentIndex;
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
        potentials = (Potential[])Arrays.addObject(potentials, newPotential);
    	return potentials.length-1;
    }

    /**
     * Removes the given potential from the list held by this type instance,
     * and returns the index that was assigned to it.  If the potential was
     * not previously added to this instance, returns -1.
     */
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
    
    private Potential[] potentials = new Potential[0];
    private int mostRecentIndex = -1;

}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.util.Arrays;

/**
 * This class stores an array of Potentials and remembers the "other" IAtomType
 * that the IPotential applies to.  PotentialMasterMonatomic holds one for each
 * IAtomType and then checks the second IAtomType held by this class to
 * identify the IPotential that applies to any given pair of IAtomTypes.
 * 
 * 1-body potentials are stored with a null IAtomType
 * 
 * @author Andrew Schultz
 */
public class PotentialArrayByType implements java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private IPotential[] potentials = new IPotential[0];
    private AtomType[] types = new AtomType[0];
    private int mostRecentIndex = -1;
    
	public PotentialArrayByType() {
		super();
	}

    /**
     * Returns the index that this type has assigned to the given
     * potential.  If potential has not been previously identified
     * to this instance, a new index will be returned.
     * @param potential The potential of interest
     * @return index associated with potential
     */
    public int getPotentialIndex(IPotential potential) {
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
    public int addPotential(IPotential newPotential, AtomType type) {
        for(mostRecentIndex=0; mostRecentIndex<potentials.length; mostRecentIndex++) {
    		if(potentials[mostRecentIndex] == newPotential) return mostRecentIndex;
    	}
        potentials = (IPotential[])Arrays.addObject(potentials, newPotential);
        types = (AtomType[]) Arrays.resizeArray(types, potentials.length);
        types[types.length-1] = type;
    	return potentials.length-1;
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

                AtomType[] newTypes = new AtomType[types.length - 1];
                System.arraycopy(types,0,newTypes,0,i);
                System.arraycopy(types,i+1,newTypes,i,types.length-i-1);
                types = newTypes;
    	    	return i;
    		}
    	}
    	return -1;
    }

    public final IPotential[] getPotentials() {
    	return potentials;
    }

    public final AtomType[] getTypes() {
        return types;
    }
}

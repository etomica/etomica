package etomica.normalmode;

import etomica.atom.Atom;

/**
 * An interface that maps real space coordinates to (per-molecule) normal mode 
 * coordinates.  For a given atom, calcU converts the molecule's real space 
 * coordinates to measures of the deviation in each of its degrees of freedom.
 * setToU places the atom into a conformation corresponding to the deviations.
 * 
 * @author Andrew Schultz
 */
public interface NormalCoordMapper {
    
    /**
     * Returns the number of normal mode coordinate dimensions needed for
     * each molecule
     */
    public int getNormalDim();
    
    /**
     * Notifies the coord wrapper how many atoms will be tracked.  The 
     * |atomCount| parameter in other methods must not exceed |numAtoms-1|.
     */
    public void setNumAtoms(int numAtoms);
    
    /**
     * Calculates the U displacements for the given atom relative to 
     * its perfect lattice value.
     * @param atom The atom of interest
     * @param nominalU The atom's nominal coordinates, from an earlier call
     * to initNominalU
     * @param u Upon return, the atom's deviation from its perfect lattice 
     * coordinates.  |u| must be of length getNormalDim()
     */
    public void calcU(Atom atom, int atomCount, double[] u);

    /**
     * Initializes the elements of nominalU to the nominal normal mode 
     * coordinates of the given atom.
     */
    public void initNominalU(Atom atom, int atomCount);

    /**
     * Set the real space coordinates for the given Atom to the normal mode
     * deltas given by u.  |u| must be of length getNormalDim()
     */
    public void setToU(Atom atom, int atomCount, double[] u);
}

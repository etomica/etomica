package etomica.lattice;

import etomica.Space;

/**
 * Interface representing a crystal basis, which is defined by a set of vectors.
 * The vectors typically represent the nominal positions of the atoms of
 * a molecule, and a molecular crystal is formed by repeating these coordinates
 * at the regular spacings defined by a Bravais lattice.
 */
 
public interface Basis {
    
	/**
	 * Number of atoms in the basis.
	 * @return int
	 */
	public int size();
    
    /**
     * Defines the nominal coordinates of the points in the basis; that
     * is, the points for a site located at the origin of the Bravais lattice.  
     * A molecular crystal would be produced by copying and translating
     * these points to the different sites of a Bravais lattice.
     */
    public Space.Vector[] positions();

}//end of Basis
package etomica.config;

import etomica.api.IBox;

/**
 * General class for assignment of molecules to positions in a box.
 * Subclasses define algorithms for the placement of molecules via their
 * implementation of the <tt>initializePositions</tt> method.  
 * <p>
 * The arrangement of atoms within the molecules is not handled by this
 * class; {@see etomica.config.Conformation}. 
 * 
 * @author David Kofke and Andrew Schultz
 */
 
public interface Configuration {

    /**
     * Defines the placement of the molecules. Atoms in all lists are 
     * considered for placement.  Subclasses define the specific algorithm
     * for placement, and for interpretation of the different lists.  
     * Some subclasses will simply place all atoms as if they are in a single
     * list, while others will place atoms in different lists in different 
     * ways, for example to form a lattice appropriate to a compound
     * (i.e., a solid-box mixture).
     * 
     * @param atomList array of list of molecules to be placed by this class
     */
    public void initializeCoordinates(IBox box);

}

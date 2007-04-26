package etomica.config;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtomGroup;
import etomica.phase.Phase;

/**
 * General class for assignment of molecules to positions in a phase.
 * Subclasses define algorithms for the placement of molecules via their
 * implementation of the <tt>initializePositions</tt> method.  
 * <p>
 * The arrangement of atoms within the molecules is not handled by this
 * class; {@see etomica.config.Conformation}. 
 * 
 * @author David Kofke and Andrew Schultz
 */
 
public abstract class Configuration implements java.io.Serializable {

    /**
     * Defines the placement of the molecules. Atoms in all lists are 
     * considered for placement.  Subclasses define the specific algorithm
     * for placement, and for interpretation of the different lists.  
     * Some subclasses will simply place all atoms as if they are in a single
     * list, while others will place atoms in different lists in different 
     * ways, for example to form a lattice appropriate to a compound
     * (i.e., a solid-phase mixture).
     * 
     * @param atomList array of list of molecules to be placed by this class
     */
    public abstract void initializeCoordinates(Phase phase);
    
    /**
     * Primary means by which this class is used to arrange molecules in the Phase.
     * Generates an array of AtomLists holding the molecules in the Phase, with
     * one list for each Species.  This array is then passed to <tt>initializePositions</tt>.
     * <p>
     * Also sets the dimensions field for this class to that of the Boundary
     * in the given Phase.
     * 
     * @param phase
     */
    protected AtomArrayList[] getMoleculeLists(Phase phase) {
        AtomArrayList speciesAgentList = phase.getSpeciesMaster().getAgentList();
        AtomArrayList[] moleculeLists = new AtomArrayList[speciesAgentList.size()];
        for (int i=0; i<speciesAgentList.size(); i++) {
            moleculeLists[i] = ((IAtomGroup)speciesAgentList.get(i)).getChildList();
        }
        return moleculeLists;
    }
    
}

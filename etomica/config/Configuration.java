package etomica.config;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.phase.Phase;
import etomica.space.Space;

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

    public Configuration(Space space) {
        this.space = space;
    }

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
    protected abstract void initializePositions(AtomArrayList[] atomList);
    
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
    public void initializeCoordinates(Phase phase) {
        setDimensions(phase.getBoundary().getDimensions().toArray());
        AtomArrayList speciesAgentList = ((AtomTreeNodeGroup)phase.getSpeciesMaster().getNode()).getChildList();
        AtomIteratorArrayListSimple speciesAgentIterator = new AtomIteratorArrayListSimple(speciesAgentList);
        AtomArrayList[] moleculeLists = new AtomArrayList[speciesAgentList.size()];
        int i=0;
        speciesAgentIterator.reset();
        while (speciesAgentIterator.hasNext()) {
            SpeciesAgent agent = (SpeciesAgent)speciesAgentIterator.nextAtom();
            moleculeLists[i++] = ((AtomTreeNodeGroup)agent.getNode()).getChildList();
        }
        initializePositions(moleculeLists);
    }
    
    /**
     * Sets the dimensions of the volume within which the molecules are arranged.
     */
    public void setDimensions(double[] dim) {
        dimensions = (double[])dim.clone();
    }

    /**
     * Returns the dimensions of the volume within which the molecules are arranged.
     */
    public double[] getDimensions() {return dimensions;}
    
    protected double[] dimensions;
    protected Space space;
    
}//end of Configuration

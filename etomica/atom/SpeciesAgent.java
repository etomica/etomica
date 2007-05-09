package etomica.atom;
import etomica.units.Dimension;
import etomica.units.Quantity;

/**
 * The SpeciesAgent is a representative of the species in each phase.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a phase.  Each phase has an agent from every species instance.
 * 
 * @author David Kofke
 */
 
public final class SpeciesAgent extends AtomGroup {

    public SpeciesAgent(AtomType type, AtomManager atomManager) {
        super(type);
        this.atomManager = atomManager;
    }

    public AtomManager getSpeciesMaster() {
        return atomManager;
    }
    
    public String signature() {
        return (atomManager != null) ? atomManager.getPhase().toString() + " " +getIndex()
                : "SpeciesAgent without phase";
    }
    
    public int getNMolecules() {return childList.size();}
            
    public IAtom addNewAtom() {
        IAtom aNew = type.getSpecies().moleculeFactory().makeAtom();
        addChildAtom(aNew);
        return aNew;
    }
    
    /**
     * Notifies this atom group that an atom has been added to it 
     * or one of its descendants.
     */
    public void addAtomNotify(IAtom childAtom) {
        atomManager.addAtomNotify(childAtom);
    }
    
    /**
     * Notifies this atom group that an atom has been removed from it or 
     * one of its descendants.
     */
    public void removeAtomNotify(IAtom childAtom) {
         atomManager.removeAtomNotify(childAtom);
    }

    /**
     * Sets the number of molecules for this species.  Makes the given number
     * of new molecules, linked-list orders and initializes them.
     * Any previously existing molecules for this species in this phase are abandoned
     * Any links to molecules of next or previous species are maintained.
     * Takes no action at all if the new number of molecules equals the existing number
     *
     * @param n  the new number of molecules for this species
     */
    public void setNMolecules(int n) {
        atomManager.notifyNewAtoms((n-getNMolecules())*type.getSpecies().moleculeFactory().getNumTreeAtoms(),
                                     (n-getNMolecules())*type.getSpecies().moleculeFactory().getNumLeafAtoms());
        if(n > childList.size()) {
            for(int i=childList.size(); i<n; i++) addNewAtom();
        }
        else if(n < childList.size()) {
            if(n < 0) {
                throw new IllegalArgumentException("Number of molecules cannot be negative");
            }
            for (int i=getChildList().size(); i>n; i--) {
                removeChildAtom(getChildList().get(i-1));
            }
        }
    }
    
    public Dimension getNMoleculesDimension() {return Quantity.DIMENSION;}

    private static final long serialVersionUID = 2L;
    private final AtomManager atomManager;
}

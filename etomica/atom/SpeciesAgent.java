package etomica.atom;
import etomica.units.Dimension;
import etomica.units.Quantity;

/**
 * The SpeciesAgent is a representative of the species in each box.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a box.  Each box has an agent from every species instance.
 * 
 * @author David Kofke
 */
 
public final class SpeciesAgent extends AtomGroup implements ISpeciesAgent {

    public SpeciesAgent(AtomType type, AtomManager atomManager) {
        super(type);
        this.atomManager = atomManager;
    }

    public AtomManager getAtomManager() {
        return atomManager;
    }
    
    public String signature() {
        return (atomManager != null) ? atomManager.getBox().toString() + " " +getIndex()
                : "SpeciesAgent without box";
    }
    
    public int getNMolecules() {return childList.getAtomCount();}
            
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

    public Dimension getNMoleculesDimension() {return Quantity.DIMENSION;}

    private static final long serialVersionUID = 2L;
    private final AtomManager atomManager;
}

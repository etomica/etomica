package etomica.atom;
import etomica.species.Species;
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

    public SpeciesAgent(AtomType type, Species species) {
        super(type);
        type.setSpecies(species);
    }
        
    public final AtomFactory moleculeFactory() {return type.getSpecies().moleculeFactory();}
      
    public int getNMolecules() {return childList.size();}
            
    public Atom addNewAtom() {
        Atom aNew = moleculeFactory().makeAtom();
        aNew.setParent(this);
        return aNew;
    }    
    
    /**
     * Overrides parent class method and terminates recursive call to identify this
     * as a constituent atom's species agent.
     */
     public final SpeciesAgent getParentSpeciesAgent() {return this;}

     /**
      * Throws a RuntimeException, because a species agent is not contained within a molecule.
      */
     public final Atom getParentMolecule() {
         throw new RuntimeException("Error:  Unexpected call to parentMolecule in SpeciesAgent");
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
        ((SpeciesMaster)parent).notifyNewAtoms((n-getNMolecules())*moleculeFactory().getNumTreeAtoms());
        if(n > childList.size()) {
            for(int i=childList.size(); i<n; i++) addNewAtom();
        }
        else if(n < childAtomCount()) {
            if(n < 0) n = 0;
            for (int i=getChildList().size(); i>n; i--) {
                getChildList().get(i-1).dispose();
            }
        }
        if (n == 0) {
            // if there are no molecules of this Species, the factory can be mutable
            // yes, this is horrible.
            SpeciesRoot speciesRoot = (SpeciesRoot)parent.getParentGroup();
            type.getSpecies().getFactory().checkMutable(speciesRoot);
        }
    }
    
    public Dimension getNMoleculesDimension() {return Quantity.DIMENSION;}

    private static final long serialVersionUID = 2L;
}

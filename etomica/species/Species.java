package etomica.species;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.IMolecule;
import etomica.simulation.SpeciesManager;

/**
 * Base implementation of ISpecies, providing basic AtomType and index
 * management.
 * @see ISpecies
 */
public abstract class Species implements java.io.Serializable, ISpecies {

    /**
     * Constructs species with molecules of the given atom type.
     */
    public Species(AtomTypeMolecule atomType) {
        this.atomType = atomType;
        atomType.setSpecies(this);
        isMutable = true;
    }
    
    public void resetIndex(SpeciesManager speciesManager) {
        ISpecies[] allSpecies = speciesManager.getSpecies();
        for (int i=0; i<allSpecies.length; i++) {
            if (allSpecies[i] == this) {
                index = i;
                return;
            }
        }
        throw new RuntimeException("I couldn't find myself.  That's bad.");
    }
    
    public int getIndex() {
        return index;
    }
    
    public abstract IMolecule makeMolecule();

    public abstract int getNumLeafAtoms();
    
    public boolean isMutable() {
        return isMutable;
    }
    
    public AtomTypeMolecule getMoleculeType() {
        return atomType;
    }
    
    /**
     * Returns a SpeciesSignature for this Species.  Subclasses must override
     * this method.
     */
    public abstract SpeciesSignature getSpeciesSignature();
    
    private static final long serialVersionUID = 2L;
    protected final AtomTypeMolecule atomType;
    protected boolean isMutable;
    protected int index;
}

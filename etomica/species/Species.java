package etomica.species;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.SpeciesManager;


 /**
  * A Species holds information about how to construct a molecule, and
  * provides for the management of molecules in a box.
  * 
  * These are the important features of a Species: <br>
  * <ol>
  * <li>It holds an AtomFactory instance that constructs molecules when
  * needed.
  * <li>It makes a SpeciesAgent class that is placed in each box (each
  * box has one species agent from each species). This agent manages the
  * molecules of that species in that box (counting, adding, removing,
  * etc.) The agent for a given box may be obtained through the getAgent
  * method. <br>
  * <li>Each Species has a unique species index assigned when it is
  * constructed. The index assignment begins at 1 and is incremented after
  * each Species construction. This index is useful when collecting things
  * in reference to the species (for example, in the use of neighbor lists).
  * </ol>
  * The number of molecules of a species in a box may be changed at run
  * time. Interactions among all molecules in a box are defined by
  * associating an intermolecular potential to one or more Species via a
  * call to the addPotential method of the PotentialMaster for the simulation.
  * 
  * @author David Kofke
  * @author C. Daniel Barnes
  * @see AtomManager
  * @see SpeciesAgent
  * @see PotentialMaster
  */
public abstract class Species implements java.io.Serializable {

    /**
     * Constructs species with molecules built by the given atom factory.
     * Species agents made by this species will have the given type for 
     * their (common) AtomType.
     */
    public Species(AtomTypeMolecule atomType) {
        this.atomType = atomType;
        atomType.setSpecies(this);
        isMutable = true;
    }
    
    public void resetIndex(SpeciesManager speciesManager) {
        Species[] allSpecies = speciesManager.getSpecies();
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
    
    /**
     * Builds and returns the atom/atomgroup made by this factory.
     * Implementation of this method in the subclass defines the 
     * product of this factory.
     */
    public abstract IMolecule makeMolecule();

    /**
     * Returns the number of leaf atoms descended from the Atom returned 
     * by makeAtom.  This will be 1 if makeAtom returns a leaf atom.
     */
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
    public SpeciesSignature getSpeciesSignature() {
        // AtomFactories can't be serialized (or rather, they'll serialize too 
        // many other things)
        return null;
    }
    
    private static final long serialVersionUID = 2L;
    protected final AtomTypeMolecule atomType;
    protected boolean isMutable;
    protected int index;
}

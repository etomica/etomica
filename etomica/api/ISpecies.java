package etomica.api;

import etomica.atom.AtomTypeMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.SpeciesManager;

/**
 * An ISpecies holds information about how to construct a molecule, and
 * creates molecules on demand.  In addition to actually creating molecules,
 * ISpecies objects are often used as keys to pass and designate which
 * species is of interest (for assigning potentials or setting the number of
 * molecules in the box, for instance).
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 * @author Andrew Schultz
 * @see IBox
 * @see PotentialMaster
 */
public interface ISpecies {

    /**
     * Informs the ISpecies that it is the child of the given SpeciesManager
     * and that it should recalculate its own index (based on its position
     * within the SpeciesManager's array of Species).
     */
    public void resetIndex(SpeciesManager speciesManager);

    /**
     * Returns the index for this ISpecies, within the context of an
     * ISimulation.  The index is the ISpecies' position in the array returned
     * by SpeciesManager.getSpecies().
     */
    public int getIndex();

    /**
     * Builds and returns the IMolecule of this ISpecies.
     */
    public IMolecule makeMolecule();

    /**
     * Returns the number of leaf atoms descended from the Atom returned 
     * by makeAtom.  This will be 1 if makeAtom returns a leaf atom.
     */
    public int getNumLeafAtoms();

    /**
     * Returns true if ISpecies methods can be called which would change the
     * definition of the Species (setFoo(foo)).  If a mutator method is called
     * when the ISpecies is not mutable, an exception is thrown.  Typically, an
     * ISpecies is mutable until a molecule is created.
     */
    public boolean isMutable();

    /**
     * Returns the atom type of molecules of this ISpecies. 
     * @return
     */
    public AtomTypeMolecule getMoleculeType();

}
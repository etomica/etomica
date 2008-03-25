package etomica.api;

import etomica.config.Conformation;
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
public interface ISpecies extends IAtomType {

    /**
     * Informs the ISpecies what its index should be.  This should only be
     * called by the SpeciesManager.
     */
    public void setIndex(int newIndex);

    /**
     * Returns the index for this ISpecies, within the context of an
     * ISimulation.  The index is the ISpecies' position in the array returned
     * by SpeciesManager.getSpecies().
     */
    public int getIndex();

    /**
     * Sets the SpeciesManager.  This is used for callbacks for notification of
     * removal and addition of child types (not that that should ever happen!)
     */
    public void setSpeciesManager(SpeciesManager newSpeciesManager);

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
     * Returns the array of child types of this group.
     */
    public IAtomTypeLeaf[] getChildTypes();

    public void addChildType(IAtomTypeLeaf newChildType);

    public void removeChildType(IAtomTypeLeaf removedType);

    /**
     * Returns the conformation used to set the standard arrangement of
     * the atoms/atom-groups produced by this factory.
     */
    public Conformation getConformation();
}
package etomica.api;

import etomica.potential.PotentialMaster;

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
     * Builds and returns the IMolecule of this ISpecies.
     */
    public IMolecule makeMolecule();

    /**
     * Returns the number of leaf atoms descended from the Atom returned 
     * by makeAtom.  This will be 1 if makeAtom returns a leaf atom.
     */
    public int getNumLeafAtoms();

    /**
     * Returns the number of child types of this group.
     */
    public int getChildTypeCount();

    /**
     * Returns the child types of this group for the specified index.
     */
    public IAtomType getChildType(int index);

    /**
     * Returns the conformation used to set the standard arrangement of
     * the atoms/atom-groups produced by this factory.
     */
    public IConformation getConformation();

    /**
     * The position definition held by the type provides an appropriate default
     * to define the position of an atom of this type. This field is set in the
     * definition of the parent species of the atom. It is null for SpeciesRoot,
     * SpeciesMaster, and SpeciesAgent atoms.
     * 
     * @return Returns the PositionDefinition for an atom of this type.
     */
    public IAtomPositionDefinition getPositionDefinition();
}

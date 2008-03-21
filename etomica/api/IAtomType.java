package etomica.api;


public interface IAtomType {

	public abstract void setIndex(int newIndex);

	public abstract int getIndex();

	/**
	 * The position definition held by the type provides an appropriate default
	 * to define the position of an atom of this type. This field is set in the
	 * definition of the parent species of the atom. It is null for SpeciesRoot,
	 * SpeciesMaster, and SpeciesAgent atoms.
	 * 
	 * @return Returns the PositionDefinition for an atom of this type.
	 */
	public abstract IAtomPositionDefinition getPositionDefinition();

	/**
	 * Sets the PositionDefinition used for this AtomType
	 */
	public abstract void setPositionDefinition(
			IAtomPositionDefinition newPositionDefinition);

	public abstract void setInteracting(boolean b);

	/**
	 * Returns true if one or more potentials are defined to act
	 * on an atom of this type.
	 */
	public abstract boolean isInteracting();

}
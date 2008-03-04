package etomica.api;

import etomica.atom.AtomPositionDefinition;

public interface IAtomType {

	public abstract void setIndex(int newIndex);

	public abstract int getIndex();

	public abstract ISpecies getSpecies();

	/**
	 * The position definition held by the type provides an appropriate default
	 * to define the position of an atom of this type. This field is set in the
	 * definition of the parent species of the atom. It is null for SpeciesRoot,
	 * SpeciesMaster, and SpeciesAgent atoms.
	 * 
	 * @return Returns the PositionDefinition for an atom of this type.
	 */
	public abstract AtomPositionDefinition getPositionDefinition();

	/**
	 * Sets the PositionDefinition used for this AtomType
	 */
	public abstract void setPositionDefinition(
			AtomPositionDefinition newPositionDefinition);

	public abstract void setInteracting(boolean b);

	/**
	 * Returns true if one or more potentials are defined to act
	 * on an atom of this type.
	 */
	public abstract boolean isInteracting();

	public abstract int compareTo(Object otherAtomType);

	public abstract String toString();

}
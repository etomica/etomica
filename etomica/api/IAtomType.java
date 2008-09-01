package etomica.api;


public interface IAtomType {

	public abstract void setIndex(int newIndex);

	public abstract int getIndex();


	public abstract void setInteracting(boolean b);

	/**
	 * Returns true if one or more potentials are defined to act
	 * on an atom of this type.
	 */
	public abstract boolean isInteracting();

}
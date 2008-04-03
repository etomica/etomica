package etomica.api;

public interface IAtomTypeSphere extends IAtomTypeLeaf {

	public abstract double getDiameter();

	/**
	 * Sets diameter of this atom and updates radius accordingly.
	 *
	 * @param d   new value for diameter
	 */
	public abstract void setDiameter(double d);

}
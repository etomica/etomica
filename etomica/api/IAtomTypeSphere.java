package etomica.api;

public interface IAtomTypeSphere extends IAtomType {

	public double getDiameter();

	/**
	 * Sets diameter of this atom and updates radius accordingly.
	 *
	 * @param d   new value for diameter
	 */
	public void setDiameter(double d);

}
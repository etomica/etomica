package etomica.api;

import etomica.atom.iterator.IteratorDirective;
import etomica.chem.models.Model;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMasterLrc;
import etomica.space.Space;

public interface IPotentialMaster {

	/**
	 * Returns the object that oversees the long-range
	 * correction zero-body potentials.
	 */
	public abstract PotentialMasterLrc lrcMaster();

	/**
	 * Returns an nBody PotentialGroup appropriate for this type of 
	 * PotentialMaster.
	 */
	public abstract PotentialGroup makePotentialGroup(int nBody);

	/**
	 * Performs the given PotentialCalculation on the atoms of the given Box.
	 * Sets the box for all molecule iterators and potentials, sets target
	 * and direction for iterators as specified by given IteratorDirective,
	 * and applies doCalculation of given PotentialCalculation with the iterators
	 * and potentials.
	 */
	public abstract void calculate(IBox box, IteratorDirective id,
			PotentialCalculation pc);

	/**
	 * Add the given Model's intramolecular potentials to this PotentialMaster
	 */
	public abstract void addModel(Model newModel);

	/**
	 * Indicates to the PotentialMaster that the given potential should apply to 
	 * the specified species.  Exception is thrown if the potential.nBody() value
	 * is different from the length of the species array.  Thus, for example, if
	 * giving a 2-body potential, then the array should contain exactly
	 * two species; the species may refer to the same instance (appropriate for an 
	 * intra-species potential, defining the iteractions between molecules of the
	 * same species).
	 */
	public abstract void addPotential(IPotential potential, ISpecies[] species);

	/**
	 * Indicates to the PotentialMaster that the given potential should apply to 
	 * the specified atom types.  The potential is assumed to be intermolecular.
	 * The given types should not include any type which is the descendent of 
	 * another.  Potential group hierarchy will be constructed as needed above
	 * the level of the given atom types.
	 * <p>
	 * The order of the elements in the atomTypes array is not relevant, and is
	 * subject to rearrangement by the method -- the array is sorted (using the compareTo
	 * method of AtomType) before doing anything else.
	 * 
	 */
	public abstract void addPotential(IPotential potential,
			IAtomType[] atomTypes);

	/**
	 * Notifies the PotentialMaster that the sub-potential has been added to 
	 * the given PotentialGroup, which is associated (but not necessarily held 
	 * by) this PotentialMaster.
	 * This method is called by PotentialGroup and should not be called in
	 * other circumstances.
	 */
	public abstract void potentialAddedNotify(IPotential subPotential,
			PotentialGroup pGroup);

	/**
	 * Returns the potential that applies to the specified types,
	 * or null of no existing potential applies.
	 */
	public abstract PotentialGroup getPotential(IAtomType[] types);

	/**
	 * Returns the AtomTypes that the given potential applies to if the given 
	 * potential is within this potential group.  If the potential is not 
	 * contained by the potential master or any PotentialGroup it holds, or 
	 * does not apply to specific AtomTypes, null is returned.
	 */
	public abstract IAtomType[] getAtomTypes(IPotential potential);

	/**
	 * Removes given potential from the group.  No error is generated if
	 * potential is not in group.
	 */
	public abstract void removePotential(IPotential potential);

	/**
	 * @return Returns enabled flag.
	 */
	public abstract boolean isEnabled();

	/**
	 * Permits enabling/disabling of all potentials.  Default is enabled (true).
	 * @param enabled flags if potentials are enabled.
	 */
	public abstract void setEnabled(boolean enabled);

	/**
	 * Indicates that the specified potential should not contribute to potential
	 * calculations. If potential is not in this group, no action is taken.
	 */
	public abstract void setEnabled(IPotential potential, boolean enabled);

	/**
	 * Returns true if the potential is in this group and has not been disabled
	 * via a previous call to setEnabled; returns false otherwise.
	 */
	public abstract boolean isEnabled(IPotential potential);

	/**
	 * @return Returns the space.
	 */
	public abstract Space getSpace();

	/**
	 * Returns an array containing all molecular Potentials.
	 */
	public abstract IPotential[] getPotentials();

}
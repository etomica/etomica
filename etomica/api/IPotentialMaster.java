package etomica.api;

import etomica.atom.iterator.IteratorDirective;
import etomica.chem.models.Model;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMasterLrc;

public interface IPotentialMaster {

	/**
	 * Returns the object that oversees the long-range
	 * correction zero-body potentials.
	 */
	public PotentialMasterLrc lrcMaster();

	/**
	 * Returns an nBody PotentialGroup appropriate for this type of 
	 * PotentialMaster.
	 */
	public PotentialGroup makePotentialGroup(int nBody);

	/**
	 * Performs the given PotentialCalculation on the atoms of the given Box.
	 * Sets the box for all molecule iterators and potentials, sets target
	 * and direction for iterators as specified by given IteratorDirective,
	 * and applies doCalculation of given PotentialCalculation with the iterators
	 * and potentials.
	 */
	public void calculate(IBox box, IteratorDirective id,
			PotentialCalculation pc);

	/**
	 * Add the given Model's intramolecular potentials to this PotentialMaster
	 */
	public void addModel(Model newModel);

	/**
	 * Indicates to the PotentialMaster that the given potential should apply to 
	 * the specified species.  Exception is thrown if the potential.nBody() value
	 * is different from the length of the species array.  Thus, for example, if
	 * giving a 2-body potential, then the array should contain exactly
	 * two species; the species may refer to the same instance (appropriate for an 
	 * intra-species potential, defining the iteractions between molecules of the
	 * same species).
	 */
	public void addPotential(IPotentialMolecular potential, ISpecies[] species);

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
	public void addPotential(IPotentialAtomic potential,
			IAtomType[] atomTypes);

	/**
	 * Notifies the PotentialMaster that the sub-potential has been added to 
	 * the given PotentialGroup, which is associated (but not necessarily held 
	 * by) this PotentialMaster.
	 * This method is called by PotentialGroup and should not be called in
	 * other circumstances.
	 */
	public void potentialAddedNotify(IPotentialAtomic subPotential,
			PotentialGroup pGroup);

	/**
	 * Returns the potential that applies to the specified types,
	 * or null of no existing potential applies.
	 */
	public PotentialGroup getPotential(ISpecies[] types);

	/**
	 * Removes given potential from the group.  No error is generated if
	 * potential is not in group.
	 */
	public void removePotential(IPotentialMolecular potential);

    /**
     * Removes given potential from the group.  No error is generated if
     * potential is not in group.
     */
    public void removePotential(IPotentialAtomic potential);

	/**
	 * @return Returns enabled flag.
	 */
	public boolean isEnabled();

	/**
	 * Permits enabling/disabling of all potentials.  Default is enabled (true).
	 * @param enabled flags if potentials are enabled.
	 */
	public void setEnabled(boolean enabled);

	/**
	 * Indicates that the specified potential should not contribute to potential
	 * calculations. If potential is not in this group, no action is taken.
	 */
	public void setEnabled(IPotentialMolecular potential, boolean enabled);

	/**
	 * Returns true if the potential is in this group and has not been disabled
	 * via a previous call to setEnabled; returns false otherwise.
	 */
	public boolean isEnabled(IPotentialMolecular potential);

    /**
     * Indicates that the specified potential should not contribute to potential
     * calculations. If potential is not in this group, no action is taken.
     */
    public void setEnabled(IPotentialAtomic potential, boolean enabled);

    /**
     * Returns true if the potential is in this group and has not been disabled
     * via a previous call to setEnabled; returns false otherwise.
     */
    public boolean isEnabled(IPotentialAtomic potential);

	/**
	 * Returns an array containing all molecular Potentials.
	 */
	public IPotentialMolecular[] getPotentials();

}
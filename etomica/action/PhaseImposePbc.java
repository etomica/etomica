/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.AtomIterator;
import etomica.AtomIteratorLeafAtoms;
import etomica.AtomIteratorMolecule;
import etomica.AtomIteratorPhaseDependent;
import etomica.Integrator;
import etomica.Phase;
import etomica.Space;

/**
 * @author kofke
 * 
 * Action that imposes the central-image effect of a phase having periodic
 * boundaries. Causes all atoms with coordinates outside the phase boundaries to
 * be moved to the central-image location (inside the boundaries).
 */

public final class PhaseImposePbc extends PhaseActionAdapter implements
		Integrator.IntervalListener {

	/**
	 * Creates the action without specifying a phase. Requires call to setPhase
	 * before action can have any effect. Default is to apply central-imaging at
	 * the atom rather than molecule level.
	 */
	public PhaseImposePbc() {
		super("Impose PBC");
		setApplyToMolecules(false);
	}
    
    public int getPriority() {return 100;}//100-199 is priority range for classes imposing PBC

	/**
	 * Creates the action ready to perform on the given phase.
	 * 
	 * @param phase
	 */
	public PhaseImposePbc(Phase phase) {
		this();
		setPhase(phase);
	}

	public void actionPerformed() {
		Space.Boundary boundary = phase.boundary();
		iterator.reset();
		while (iterator.hasNext()) {
			boundary.centralImage(iterator.nextAtom().coord);
		}
	}

	public void intervalAction(Integrator.IntervalEvent evt) {
		actionPerformed();
	}

	public void setPhase(Phase phase) {
		super.setPhase(phase);
		iterator.setPhase(phase);
	}

	/**
	 * Returns the iterator that gives the atoms to which central imaging is
	 * applied.
	 * 
	 * @return AtomIteratorList
	 */
	public AtomIterator getIterator() {
		return iterator;
	}

	/**
	 * Sets the iterator the gives the atoms for central imaging. Normally this
	 * does not need to be set, but if central imaging scheme is desired for
	 * application at a level other than the molecules or the leaf atoms, or if
	 * it is to be applied only to a subet of the atoms in a phase, this can be
	 * invoked by setting this iterator appropriately.
	 * 
	 * @param iterator
	 *            The iterator to set
	 */
	public void setIterator(AtomIteratorPhaseDependent iterator) {
		this.iterator = iterator;
		iterator.setPhase(phase);
	}

	/**
	 * Returns the value of applyToMolecules.
	 * 
	 * @return boolean
	 */
	public boolean isApplyToMolecules() {
		return applyToMolecules;
	}

	/**
	 * Sets a flag indicating whether periodic boundaries are applied to the
	 * molecules (true), or to the atoms (false). If applied to the atoms (the
	 * default case), then central imaging is done to each atom individually,
	 * which could cause a molecule to be split, with some of its atoms on one
	 * edge of the simulation box, and others on the other edge. If applied to
	 * molecules, the entire molecule will be shifted as a whole when enforcing
	 * central imaging.
	 * 
	 * @param applyToMolecules
	 *            The new value of the flag.
	 */
	public void setApplyToMolecules(boolean applyToMolecules) {
		this.applyToMolecules = applyToMolecules;
		if (applyToMolecules)
			iterator = new AtomIteratorMolecule(phase);
		else
			iterator = new AtomIteratorLeafAtoms(phase);
	}

	private AtomIteratorPhaseDependent iterator;

	private boolean applyToMolecules;

}
/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Integrator;
import etomica.Phase;
import etomica.Space;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceCOM;
import etomica.space.Boundary;
import etomica.space.Vector;

/**
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
		Boundary boundary = phase.boundary();
		iterator.reset();
        if (applyToMolecules) {
            while (iterator.hasNext()) {
                Atom molecule = iterator.nextAtom();
                Vector shift = boundary.centralImage(moleculePosition.position(molecule));
                if (!shift.isZero()) {
                    translator.setTranslationVector(shift);
                    moleculeTranslator.actionPerformed(molecule);
                }
            }
        }
        else {
            while (iterator.hasNext()) {
                Atom atom = iterator.nextAtom();
                Vector shift = boundary.centralImage(atom.coord.position());
                if (!shift.isZero()) {
                    atom.coord.position().PE(shift);
                }
            }
        }
	}

	public void intervalAction(Integrator.IntervalEvent evt) {
		actionPerformed();
	}

	public void setPhase(Phase phase) {
		super.setPhase(phase);
		iterator.setPhase(phase);
        if (space != phase.space()) {
            space = phase.space();
            translator = new AtomActionTranslateBy(phase.space());
            moleculeTranslator = new AtomGroupAction(translator);
            //XXX shouldn't clobber user-set moleculePosition, but DataSourceCOM needs the space
            moleculePosition = new DataSourceCOM(phase.space());
        }
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
			iterator = new AtomIteratorAllMolecules(phase);
		else
			iterator = new AtomIteratorLeafAtoms(phase);
	}

    public void setMoleculePositionDefintion(AtomPositionDefinition positionDefinition) {
        moleculePosition = positionDefinition;
    }
    public AtomPositionDefinition getMoleculePositionDefintion() {
        return moleculePosition;
    }
    
	private AtomIteratorPhaseDependent iterator;
    private AtomActionTranslateBy translator;
    private AtomGroupAction moleculeTranslator;
    private AtomPositionDefinition moleculePosition;
    private Space space;

	private boolean applyToMolecules;

}
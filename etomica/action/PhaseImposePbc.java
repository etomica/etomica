package etomica.action;

import etomica.atom.AtomPositionCOM;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Action that imposes the central-image effect of a phase having periodic
 * boundaries. Causes all atoms with coordinates outside the phase boundaries to
 * be moved to the central-image location (inside the boundaries).
 */

public final class PhaseImposePbc extends PhaseActionAdapter {

	/**
	 * Creates the action without specifying a phase. Requires call to setPhase
	 * before action can have any effect. Default is to apply central-imaging at
	 * the atom rather than molecule level.
	 */
	public PhaseImposePbc() {
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
		Boundary boundary = phase.getBoundary();
		iterator.reset();
        if (applyToMolecules) {
            for (IAtom molecule = iterator.nextAtom(); molecule != null;
                 molecule = iterator.nextAtom()) {
                IVector shift = boundary.centralImage(moleculePosition.position(molecule));
                if (!shift.isZero()) {
                    translator.setTranslationVector(shift);
                    moleculeTranslator.actionPerformed(molecule);
                }
            }
        }
        else {
            for (IAtomPositioned atom = (IAtomPositioned)iterator.nextAtom(); atom != null;
                 atom = (IAtomPositioned)iterator.nextAtom()) {
                IVector shift = boundary.centralImage(atom.getPosition());
                if (!shift.isZero()) {
                    atom.getPosition().PE(shift);
                }
            }
        }
	}

	public void setPhase(Phase phase) {
		super.setPhase(phase);
		iterator.setPhase(phase);
        if (space != phase.getSpace()) {
            space = phase.getSpace();
            translator = new AtomActionTranslateBy(phase.getSpace());
            moleculeTranslator = new AtomGroupAction(translator);
            if (moleculePosition == null) {
                moleculePosition = new AtomPositionCOM(phase.getSpace());
            }
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
			iterator = new AtomIteratorAllMolecules();
		else
			iterator = new AtomIteratorLeafAtoms();
        if (phase != null) {
            iterator.setPhase(phase);
        }
	}

    public void setMoleculePositionDefintion(AtomPositionDefinition positionDefinition) {
        moleculePosition = positionDefinition;
    }
    public AtomPositionDefinition getMoleculePositionDefintion() {
        return moleculePosition;
    }
    
    private static final long serialVersionUID = 1L;
	private AtomIteratorPhaseDependent iterator;
    private AtomActionTranslateBy translator;
    private AtomGroupAction moleculeTranslator;
    private AtomPositionDefinition moleculePosition;
    private Space space;

	private boolean applyToMolecules;
}
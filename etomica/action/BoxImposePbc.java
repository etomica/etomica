package etomica.action;

import etomica.api.IVector;
import etomica.atom.AtomPositionCOM;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.Space;

/**
 * Action that imposes the central-image effect of a box having periodic
 * boundaries. Causes all atoms with coordinates outside the box boundaries to
 * be moved to the central-image location (inside the boundaries).
 */

public class BoxImposePbc extends BoxActionAdapter {

	/**
	 * Creates the action without specifying a box. Requires call to setBox
	 * before action can have any effect. Default is to apply central-imaging at
	 * the atom rather than molecule level.
	 */
	public BoxImposePbc(Space space) {
		setApplyToMolecules(false);
		this.space = space;
	}
    
    public int getPriority() {return 100;}//100-199 is priority range for classes imposing PBC

	/**
	 * Creates the action ready to perform on the given box.
	 * 
	 * @param box
	 */
	public BoxImposePbc(Box box, Space space) {
		this(space);
		setBox(box);
	}

	public void actionPerformed() {
		Boundary boundary = box.getBoundary();
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

	public void setBox(Box box) {
		super.setBox(box);
		iterator.setBox(box);
        if (space.D() != box.getBoundary().getDimensions().getD()) {
            throw new IllegalArgumentException("Cannot change dimension of BoxImosePbc");
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
	 * it is to be applied only to a subet of the atoms in a box, this can be
	 * invoked by setting this iterator appropriately.
	 * 
	 * @param iterator
	 *            The iterator to set
	 */
	public void setIterator(AtomIteratorBoxDependent iterator) {
		this.iterator = iterator;
		iterator.setBox(box);
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
        if (box != null) {
            iterator.setBox(box);
        }
	}

    public void setMoleculePositionDefintion(AtomPositionDefinition positionDefinition) {
        moleculePosition = positionDefinition;
    }
    public AtomPositionDefinition getMoleculePositionDefintion() {
        return moleculePosition;
    }
    
    private static final long serialVersionUID = 1L;
	private AtomIteratorBoxDependent iterator;
    private AtomActionTranslateBy translator;
    private AtomGroupAction moleculeTranslator;
    private AtomPositionDefinition moleculePosition;
    private Space space;

	private boolean applyToMolecules;
}
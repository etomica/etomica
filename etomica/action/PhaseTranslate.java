/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.Space;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.space.Vector;

/**
 * Action that moves (translates) all the atoms in a phase by a specified
 * (vector) amount.
 */

public final class PhaseTranslate extends PhaseActionAdapter implements Undoable {

	private final Vector translationVector;
    private final AtomIteratorLeafAtoms iterator;

	/**
	 * Constructor requires space instance to define a translation vector.
	 * 
	 * @param space
	 */
	public PhaseTranslate(Space space) {
		super("PhaseTranslate");
		translationVector = space.makeVector();
        iterator = new AtomIteratorLeafAtoms();
        
	}

	/**
	 * Sets vector specifying how atoms are to be moved. Vector is copied so to
	 * update it is necessary to invoke this method (cannot effect change simply
	 * by changing the coordinates of vector previously passed to this method).
	 */
	public void setTranslationVector(Vector v) {
		translationVector.E(v);
	}

	/**
	 * @return the translation vector. Changes to coordinates of the returned
	 *         vector will change the behavior of this action.
	 */
	public Vector getTranslationVector() {
		return translationVector;
	}

	/**
	 * Translates all atoms by the amount of the current tranlationVector.
	 */
	public void actionPerformed() {
		if (phase == null)
			return;
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {
            iterator.nextAtom().coord.position().PE(translationVector);
        }
	}

	/**
	 * Reverses the outcome of actionPerformed, moving the atoms an amount
	 * opposite to translationVector. Assumes translationVector and atom
	 * positions have not been changed since actionPerformed.
	 */
	public void undo() {
		translationVector.TE(-1);
		actionPerformed();
		translationVector.TE(-1);
	}
}
/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.AtomIteratorAllMolecules;
import etomica.Phase;

/**
 * Performs actions that cause volume of system to expand, with molecule
 * positions scaled to keep them in the same relative positions. Inflation is
 * done isotropically, equal in all coordinate directions.
 */
public final class PhaseInflate extends PhaseActionAdapter implements Undoable {

	/**
	 * Constructs action with a default scale of 1.0.  Requires call
	 * to setPhase before action can have any effect.
	 */
	public PhaseInflate() {
		super("Inflate");
		moleculeIterator = new AtomIteratorAllMolecules();
		setScale(1.0);
	}

	/**
	 * Constructs action ready to be performed on the given phase. 
	 */
	public PhaseInflate(Phase phase) {
		this();
		setPhase(phase);
	}

	/**
	 * Sets the scale defining the amount of inflation. A value of 1.0 causes no
	 * change, while a value greater than 1.0 expands the phase, and a value
	 * less than 1.0 contracts the phase. Coordinates and boundary dimensions
	 * are all multiplied by the scale when action is performed. A zero or
	 * negative scale throws an IllegalArgumentException.
	 */
	public void setScale(double scale) {
		if (scale <= 0.0)
			throw new IllegalArgumentException(
					"Cannot have zero or negative scaling in PhaseInflate");
		this.scale = scale;
	}

	/**
	 * @return Current value of the inflation scale.
	 */
	public double getScale() {
		return scale;
	}

	/**
	 * Sets the phase to which the action will be applied.
	 */
	public void setPhase(Phase phase) {
		super.setPhase(phase);
		moleculeIterator.setPhase(phase);
	}

	/**
	 * Performs isotropic inflation.
	 */
	public void actionPerformed() {
		if(phase == null) return;
		phase.boundary().inflate(scale);
		moleculeIterator.reset();
		while (moleculeIterator.hasNext()) {
			moleculeIterator.nextAtom().coord.inflate(scale);
		}
	}

	/**
	 * Reverses the action of the inflation by performing the
	 * action with a scale given the by the reciprocal of the 
	 * current scale.  Value of scale is not changed as a result.
	 */
	public void undo() {
		setScale(1.0 / scale);
		actionPerformed();
		setScale(1.0 / scale);
	}

	private AtomIteratorAllMolecules moleculeIterator;

	private double scale = 1.0;

}
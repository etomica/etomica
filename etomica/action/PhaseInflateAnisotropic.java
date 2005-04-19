/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.Phase;
import etomica.PhaseEvent;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.space.Vector;

/**
 * Performs actions that cause volume of system to expand, with particle
 * positions scaled to keep them in the same relative positions. Inflation is
 * done anisotropically, so that each dimension can be scaled differently.
 */
public class PhaseInflateAnisotropic extends PhaseActionAdapter implements
		Undoable {

	/**
	 * Constructs action ready to be performed on the given phase.
	 */
	public PhaseInflateAnisotropic(Phase phase) {
		this(phase.space());
		setPhase(phase);
        inflateEvent = new PhaseEvent(this,PhaseEvent.BOUNDARY_INFLATE);
	}

	/**
	 * Constructs action with a default scale of 1.0 in each direction. Requires
	 * call to setPhase before action can have any effect.
	 */
	public PhaseInflateAnisotropic(Space space) {
		super("Anisotropic Inflate");
		scale = space.makeVector();
		temp = space.makeVector();
		oldDimensions = space.makeVector();
	}

	/**
	 * Sets the scale defining the amount of inflation. A value of 1.0 for any
	 * element causes no change in scale in the corresponding spatial dimension,
	 * while a value greater than 1.0 expands the phase in the corresponding
	 * direction, and a value less than 1.0 contracts. Coordinates and boundary
	 * dimensions are all multiplied by the scale when action is performed. A
	 * zero or negative scale for any element throws an
	 * IllegalArgumentException.
	 */
	public void setScale(Vector scale) {
		for (int i = 0; i < scale.D(); i++) {
			if (scale.x(i) <= 0.0)
				throw new IllegalArgumentException(
						"Cannot have zero or negative scaling in PhaseInflateAnisotropic");
		}
		this.scale.E(scale);
	}

	/**
	 * @return a clone of the current scale vector.
	 */
	public Vector getScale() {
		return (Vector) scale.clone();
	}

	/**
	 * Sets the phase to which the action will be applied.
	 */
	public void setPhase(Phase phase) {
		super.setPhase(phase);
		moleculeIterator.setPhase(phase);
	}

	/**
	 * Performs anisotropic inflation using current phase and scale.
	 */
	public void actionPerformed() {
        Vector dimensions = phase.boundary().dimensions();
        dimensions.TE(scale);
        phase.setDimensions(dimensions);
        phase.boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
		moleculeIterator.reset();
		while (moleculeIterator.hasNext()) {
			moleculeIterator.nextAtom().coord.position().TE(scale);
		}
	}

	/**
	 * Reverses the action of the inflation by performing the action with a
	 * scale given the by the reciprocal of the current scale. Value of scale is
	 * not changed as a result.
	 */
	public void undo() {
		temp.E(1.0);
		temp.DE(scale);
		scale.E(temp);
		actionPerformed();
		temp.E(1.0);
		temp.DE(scale);
		scale.E(temp);
	}

	private AtomIteratorAllMolecules moleculeIterator;
	protected Vector scale;
	protected transient Vector temp;
	private transient Vector oldDimensions;
    private PhaseEvent inflateEvent;

}
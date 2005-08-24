/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.Phase;
import etomica.atom.Atom;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Performs actions that cause volume of system to expand, with molecule
 * positions scaled to keep them in the same relative positions. Inflation can
 * be isotropically or anisotropically.
 */
public final class PhaseInflate extends PhaseActionAdapter implements Undoable {

    /**
     * Constructs action with a default scale of 1.0.  Requires call
     * to setPhase before action can have any effect.
     */
    public PhaseInflate(Space space) {
        super("Inflate");
        translator = new AtomActionTranslateBy(space);
        groupScaler = new AtomGroupAction(translator);
        moleculeIterator = new AtomIteratorAllMolecules();
        moleculeCenter = new AtomPositionGeometricCenter(space);
        scaleVector = space.makeVector();
        setScale(1.0);
    }

    /**
     * Constructs action ready to be performed on the given phase. 
     */
    public PhaseInflate(Phase phase) {
        this(phase.space());
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
        scaleVector.E(scale);
    }

    /**
     * @return Current value of the inflation scale.
     * Assumes isotropic scaling had been set eariler.
     */
    public double getScale() {
        return scaleVector.x(0);
    }

    /**
     * Sets the scale defining the amount of inflation for each dimension. 
     * A value of 1.0 causes no change, while a value greater than 1.0 expands 
     * the phase, and a value less than 1.0 contracts the phase. Coordinates 
     * and boundary dimensions are all multiplied by the scale when action is 
     * performed. A zero or negative scale throws an IllegalArgumentException.
     */
    public void setVectorScale(Vector scale) {
        if (scale.min() <= 0.0) throw new IllegalArgumentException(
                   "Cannot have zero or negative scaling in PhaseInflate");
        scaleVector.E(scale);
    }
    
    /**
     * Returns the current value of the inflation scale in each dimension.
     */
    public Vector getVectorScale() {
        return scaleVector;
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
        Vector dimensions = phase.boundary().dimensions();
        dimensions.TE(scaleVector);
        phase.setDimensions(dimensions);
        moleculeIterator.reset();
        // substract 1 from each dimension so that multiplying by it yields
        // the amount each coordinate is to be translated *by* (not to).
        scaleVector.PE(-1);
        Vector translationVector = translator.getTranslationVector();
        while (moleculeIterator.hasNext()) {
            Atom molecule = moleculeIterator.nextAtom();
            translationVector.E(moleculeCenter.position(molecule));
            translationVector.TE(scaleVector);
            groupScaler.actionPerformed(molecule);
        }
        // undo previous subtraction
        scaleVector.PE(1);
    }

    /**
     * Reverses the action of the inflation by performing the
     * action with a scale given the by the reciprocal of the 
     * current scale.  Value of scale is not changed as a result.
     */
    public void undo() {
        for (int i=0; i<scaleVector.D(); i++) {
            scaleVector.setX(i,1/scaleVector.x(i));
        }
        actionPerformed();
        for (int i=0; i<scaleVector.D(); i++) {
            scaleVector.setX(i,1/scaleVector.x(i));
        }
    }

    private final AtomIteratorAllMolecules moleculeIterator;
    private final AtomActionTranslateBy translator;
    private final AtomGroupAction groupScaler;
    private final Vector scaleVector;
    private final AtomPositionGeometricCenter moleculeCenter;
}
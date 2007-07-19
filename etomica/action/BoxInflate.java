/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.box.Box;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.util.Debug;

/**
 * Performs actions that cause volume of system to expand, with molecule
 * positions scaled to keep them in the same relative positions. Inflation can
 * be isotropically or anisotropically.
 */
public class BoxInflate extends BoxActionAdapter implements Undoable {

    /**
     * Constructs action with a default scale of 1.0.  Requires call
     * to setBox before action can have any effect.
     */
    public BoxInflate(Space space) {
        translator = new AtomActionTranslateBy(space);
        groupScaler = new AtomGroupAction(translator);
        moleculeIterator = new AtomIteratorAllMolecules();
        moleculeCenter = new AtomPositionGeometricCenter(space);
        scaleVector = space.makeVector();
        setScale(1.0);
    }

    /**
     * Constructs action ready to be performed on the given box. 
     */
    public BoxInflate(Box box) {
        this(box.getSpace());
        setBox(box);
    }

    /**
     * Sets the scale defining the amount of inflation. A value of 1.0 causes no
     * change, while a value greater than 1.0 expands the box, and a value
     * less than 1.0 contracts the box. Coordinates and boundary dimensions
     * are all multiplied by the scale when action is performed. A zero or
     * negative scale throws an IllegalArgumentException.
     */
    public void setScale(double scale) {
        if (scale <= 0.0)
            throw new IllegalArgumentException(
                    "Cannot have zero or negative scaling in BoxInflate");
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
     * the box, and a value less than 1.0 contracts the box. Coordinates 
     * and boundary dimensions are all multiplied by the scale when action is 
     * performed. A zero or negative scale throws an IllegalArgumentException.
     */
    public void setVectorScale(IVector scale) {
        scaleVector.E(scale);
    }
    
    /**
     * Returns the current value of the inflation scale in each dimension.
     */
    public IVector getVectorScale() {
        return scaleVector;
    }
    
    /**
     * Sets the box to which the action will be applied.
     */
    public void setBox(Box box) {
        super.setBox(box);
        moleculeIterator.setBox(box);
    }

    /**
     * Performs isotropic inflation.
     */
    public void actionPerformed() {
        if(box == null) return;
        IVector dimensions = box.getBoundary().getDimensions();
        dimensions.TE(scaleVector);
        box.getBoundary().setDimensions(dimensions);
        moleculeIterator.reset();
        // substract 1 from each dimension so that multiplying by it yields
        // the amount each coordinate is to be translated *by* (not to).
        scaleVector.PE(-1);
        IVector translationVector = translator.getTranslationVector();
        
        for (IAtom molecule = moleculeIterator.nextAtom(); molecule != null;
             molecule = moleculeIterator.nextAtom()) {
            translationVector.E(moleculeCenter.position(molecule));
            translationVector.TE(scaleVector);
            groupScaler.actionPerformed(molecule);
        }
        
        // undo previous subtraction
        scaleVector.PE(1);
        
        //XXX pretend we're just setting the dimensions now.  the only effect
        // here is to fire an event that notifies things (NeighborCellManager)
        // that the box changed and it should (perhaps) resize the lattice and
        // reassign atoms to cells
        box.setDimensions(dimensions);
    }

    /**
     * Reverses the action of the inflation by performing the
     * action with a scale given the by the reciprocal of the 
     * current scale.  Value of scale is not changed as a result.
     */
    public void undo() {
        for (int i=0; i<scaleVector.getD(); i++) {
            scaleVector.setX(i,1/scaleVector.x(i));
        }
        actionPerformed();
        for (int i=0; i<scaleVector.getD(); i++) {
            scaleVector.setX(i,1/scaleVector.x(i));
        }
    }

    private static final long serialVersionUID = 1L;
    protected final AtomIteratorAllMolecules moleculeIterator;
    protected final AtomActionTranslateBy translator;
    protected final AtomGroupAction groupScaler;
    protected final IVector scaleVector;
    protected final AtomPositionGeometricCenter moleculeCenter;
    
}
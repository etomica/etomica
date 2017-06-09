/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.MoleculePositionGeometricCenter;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * Performs actions that cause volume of system to expand or contract, with
 * molecule positions scaled to keep them in the same relative positions.
 * Inflation can be isotropically or anisotropically.
 */
public class BoxInflate extends BoxActionAdapter implements Undoable {

    /**
     * Constructs action with a default scale of 1.0.  Requires call
     * to setBox before action can have any effect.
     */
    public BoxInflate(Space space) {

        translator = new AtomActionTranslateBy(space);
        groupScaler = new MoleculeChildAtomAction(translator);
        moleculeCenter = new MoleculePositionGeometricCenter(space);
        scaleVector = space.makeVector();
        dimVector = space.makeVector();
        setScale(1.0);
    }

    /**
     * Constructs action ready to be performed on the given box. 
     */
    public BoxInflate(Box box, Space space) {
        this(space);
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
        return scaleVector.getX(0);
    }

    /**
     * Sets the scale defining the amount of inflation for each dimension. 
     * A value of 1.0 causes no change, while a value greater than 1.0 expands 
     * the box, and a value less than 1.0 contracts the box. Coordinates 
     * and boundary dimensions are all multiplied by the scale when action is 
     * performed. A zero or negative scale throws an IllegalArgumentException.
     */
    public void setVectorScale(Vector scale) {
        scaleVector.E(scale);
    }
    
    /**
     * Returns the current value of the inflation scale in each dimension.
     */
    public Vector getVectorScale() {
        return scaleVector;
    }
    
    /**
     * Sets the action to change the box dimensions isotropically to achieve
     * the given density.
     */
    public void setTargetDensity(double newTargetDensity) {
        double vNew = box.getMoleculeList().getMoleculeCount()/newTargetDensity;
        double scale = Math.pow(vNew/box.getBoundary().volume(), 1.0/scaleVector.getD());
        setScale(scale);
    }
    
    /**
     * Returns the target density of the action.
     */
    public double getTargetDensity() {
        double rho = box.getMoleculeList().getMoleculeCount()/box.getBoundary().volume();
        for (int i=0; i<scaleVector.getD(); i++) {
            rho *= scaleVector.getX(i);
        }
        return rho;
    }

    /**
     * Performs boundary dimension change
     */
    public void actionPerformed() {
        if(box == null) return;
        // substract 1 from each dimension so that multiplying by it yields
        // the amount each coordinate is to be translated *by* (not to).
        scaleVector.PE(-1);
        Vector translationVector = translator.getTranslationVector();

        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            translationVector.E(moleculeCenter.position(molecule));
            translationVector.TE(scaleVector);
            groupScaler.actionPerformed(molecule);
        }
        
        // undo previous subtraction
        scaleVector.PE(1);

        // actually change the boundary.  With cell/neighbor listing, this
        // triggers cell reassignment, which would fail if we shrink before
        // shifting atoms.  If we grow first, cell assignment would succeed,
        // but we'd need to reassign again after scaling the atom positions.
        dimVector.E(box.getBoundary().getBoxSize());
        dimVector.TE(scaleVector);
        box.getBoundary().setBoxSize(dimVector);
    }

    /**
     * Reverses the action of the inflation by performing the
     * action with a scale given the by the reciprocal of the 
     * current scale.  Value of scale is not changed as a result.
     */
    public void undo() {
        for (int i=0; i<scaleVector.getD(); i++) {
            scaleVector.setX(i,1/scaleVector.getX(i));
        }
        actionPerformed();
        for (int i=0; i<scaleVector.getD(); i++) {
            scaleVector.setX(i,1/scaleVector.getX(i));
        }
    }

    private static final long serialVersionUID = 1L;
    protected final AtomActionTranslateBy translator;
    protected final MoleculeChildAtomAction groupScaler;
    protected final Vector scaleVector, dimVector;
    protected final MoleculePositionGeometricCenter moleculeCenter;
    
}

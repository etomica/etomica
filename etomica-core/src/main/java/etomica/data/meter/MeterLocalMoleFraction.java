/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Fraction;

/**
 * Meter for measurement of species mole fraction within a specified subvolume
 */
public class MeterLocalMoleFraction extends DataSourceScalar {

    public MeterLocalMoleFraction(Space space, Box _box) {
        super("Local Mole Fraction",Fraction.DIMENSION);
        if(!(_box.getBoundary() instanceof Boundary)) {
        	throw new RuntimeException("The box boundary must be a subclass of etomica.Space.Boundary");
        }
        box = _box;
        tempVec = space.makeVector();
        shapeOrigin = space.makeVector();
        iterator.setBox(box);
        setShape(((Boundary)box.getBoundary()).getShape());
    }

    /**
     * Sets the subvolume shape for the mole fraction calculation.
     */
    public void setShape(Polytope shape) {
        this.shape = shape;
    }

    /**
     * Returns the subvolume shape for the mole fraction calculation.
     */
    public Polytope getShape() {
        return shape;
    }
    
    /**
     * Sets the origin of the subvolume (Polytopes typically have their center
     * at 0).  The shape origin is 0 by default.
     */
    public void setShapeOrigin(Vector newShapeOrigin) {
        shapeOrigin = newShapeOrigin;
    }
    
    /**
     * Returns the origin of the subvolume.
     */
    public Vector getShapeOrigin() {
        return shapeOrigin;
    }
    
    /**
     * Accessor method to set which species mole-fraction or molar-density is averaged
     * To set to total number density, invoke with static ALL_SPECIES field as argument
     */
    public final void setSpecies(ISpecies s) {species = s;}
    /**
     * Accessor method to get the current value of the species index
     *
     */
    public final ISpecies getSpecies() {return species;}
    
    /**
     * @return the current value of the local density or local mole fraction
     */
    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        int totalSum = 0, speciesSum = 0;
        iterator.reset();
        for (IAtom a = iterator.nextAtom(); a != null;
             a = iterator.nextAtom()) {
            tempVec.Ev1Mv2(a.getPosition(), shapeOrigin);
            if(shape.contains(tempVec)) {
                totalSum++;
                if(a.getType().getSpecies() == species) speciesSum++;
            }
        }
        if(totalSum == 0) return Double.NaN;
        return (double)speciesSum/(double)totalSum;
    }
    
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
/*
    public void setBox(Box newBox) {
        box = newBox;
        tempVec = space.makeVector();
        shapeOrigin = space.makeVector();
        iterator.setBox(box);
        if (shape == null) {
            setShape(box.getBoundary().getShape());
        }
    }
*/
    /**
     * @return Returns the iterator.
     */
    public AtomIteratorBoxDependent getIterator() {
        return iterator;
    }
    /**
     * @param iterator The iterator to set.
     */
    public void setIterator(AtomIteratorBoxDependent iterator) {
        this.iterator = iterator;
    }

    private static final long serialVersionUID = 1L;
    private Box box;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private ISpecies species;
    private AtomIteratorBoxDependent iterator = new AtomIteratorLeafAtoms();
    private Polytope shape;
    private Vector shapeOrigin;
    private Vector tempVec;
}

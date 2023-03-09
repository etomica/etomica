/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Fraction;

/**
 * Meter for measurement of species mole fraction within a specified subvolume
 */
public class MeterLocalMoleFraction extends DataSourceScalar {

    public MeterLocalMoleFraction(Space space, Box _box) {
        super("Local Mole Fraction",Fraction.DIMENSION);
        box = _box;
        tempVec = space.makeVector();
        shapeOrigin = space.makeVector();
        setShape(box.getBoundary().getShape());
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
        for (IAtom a : box.getLeafList()) {
            tempVec.Ev1Mv2(a.getPosition(), shapeOrigin);
            if(shape.contains(tempVec)) {
                totalSum++;
                if(a.getType().getSpecies() == species) speciesSum++;
            }
        }
        if(totalSum == 0) return Double.NaN;
        return (double)speciesSum/(double)totalSum;
    }
    
    private final Box box;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private ISpecies species;
    private Polytope shape;
    private Vector shapeOrigin;
    private final Vector tempVec;
}

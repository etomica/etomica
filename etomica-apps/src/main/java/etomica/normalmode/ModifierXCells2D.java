/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.modifier.Modifier;
import etomica.species.ISpecies;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

/**
 * Modifier class that enables change of the number of cells in 2D
 * 
 */
public class ModifierXCells2D implements Modifier, java.io.Serializable {

    /**
     * 
     */
    public ModifierXCells2D(Box box, ISpecies species, int y) {
        this.box = box;
        this.species = species;
        this.yCell = y;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        previousValue = mostRecentValue;
        mostRecentValue = (int)d;
        box.setNMolecules(species, 2*(int)d*yCell);
    }

    public double getValue() {
        return box.getNMolecules(species)/(2*yCell);
    }

    public Dimension getDimension() {
        return Null.DIMENSION;
    }
    
    public String getLabel() {
        return "x-Cell Number";
    }
    
 
    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final ISpecies species;
    protected int mostRecentValue, previousValue, yCell;
}

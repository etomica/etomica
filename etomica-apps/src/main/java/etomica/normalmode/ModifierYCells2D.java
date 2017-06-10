/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.modifier.Modifier;
import etomica.species.ISpecies;
import etomica.units.Dimension;
import etomica.units.Null;

/**
 * Modifier class that enables change of the number of cells in 2D
 * 
 */
public class ModifierYCells2D implements Modifier, java.io.Serializable {

    /**
     * 
     */
    public ModifierYCells2D(Box box, ISpecies species, int x) {
        this.box = box;
        this.species = species;
        this.xCell = x;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        previousValue = mostRecentValue;
        mostRecentValue = (int)d;
        
        box.setNMolecules(species, 2*(int)d*xCell);
    }

    public double getValue() {
        return box.getNMolecules(species)/(2*xCell);
    }

    public Dimension getDimension() {
        return Null.DIMENSION;
    }
    
    public String getLabel() {
        return "y-Cell Number";
    }
    
 
    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final ISpecies species;
    protected int mostRecentValue, previousValue, xCell;
}

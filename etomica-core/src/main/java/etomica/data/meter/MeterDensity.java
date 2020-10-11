/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;

/**
 * Meter for measurement of the total molecule number density in a box
 * Molecule number density is defined (number of molecules)/(volume of box)
 */
public class MeterDensity extends DataSourceScalar {
    
    public MeterDensity(Box box) {
        super("Number Density",new DimensionRatio(Quantity.DIMENSION, Volume.dimension(box.getSpace().D())));
        this.box = box;
    }

    public void setSpecies(ISpecies s) {
        species = s;
    }
    public ISpecies getSpecies() {
    	return species;
    }

    public double getDataAsScalar() {
        return (species == null ?
        			box.getMoleculeList().size() :
        			box.getNMolecules(species))
				/box.getBoundary().volume();
    }
    
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    private final Box box;
    private ISpecies species;
}

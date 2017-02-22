/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modifier;

import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.units.Dimension;
import etomica.units.Quantity;

/**
 * Modifier class that enables change of the number of molecules of a particular species
 * in a particular box.
 */
public class ModifierNMolecule implements Modifier, java.io.Serializable {

    /**
     * @param speciesAgent Agent of the affected species in the affected box.
     * Cannot be changed after construction.
     */
    public ModifierNMolecule(IBox box, ISpecies species) {
        this.box = box;
        this.species = species;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        box.setNMolecules(species, (int) d);
    }

    public double getValue() {
        return box.getNMolecules(species);
    }

    public Dimension getDimension() {
        return Quantity.DIMENSION;
    }
    
    public String getLabel() {
        return species + " molecules";
    }
    
    private static final long serialVersionUID = 1L;
    protected final IBox box;
    protected final ISpecies species;
}

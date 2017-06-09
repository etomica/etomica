/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.api.ISpecies;
import etomica.data.DataSourceMolecular;
import etomica.data.DataSourceScalar;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.units.Quantity;

/**
 * Meter for recording the total number of molecules in the box
 */
public class MeterNMolecules extends DataSourceScalar implements DataSourceMolecular {
    
    private static final long serialVersionUID = 1L;
    private ISpecies species;
    
    public MeterNMolecules() {
        super("Molecules",Quantity.DIMENSION);
    }

    public void setSpecies(ISpecies s) {species = s;}
    public ISpecies getSpecies() {return species;}

    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        return (species == null) ? box.getMoleculeList().getMoleculeCount(): box.getNMolecules(species);
    }
    
    public IData getData(IMolecule atom) {
        data.x = (species == null || (atom.getType() == species)) ? 1 : 0;
        return data;
    }
    
    public IEtomicaDataInfo getMoleculeDataInfo() {
        return dataInfo;
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
    public void setBox(Box box) {
        this.box = box;
    }

    private Box box;
}

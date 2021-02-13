/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceMolecular;
import etomica.data.DataSourceScalar;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.molecule.IMolecule;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Temperature;

/**
 * Meter for recording the total number of molecules in the box
 */
public class MeterMoleculeTemperature extends DataSourceScalar implements DataSourceMolecular {

    private ISpecies species;

    public MeterMoleculeTemperature() {
        super("Temperature", Temperature.DIMENSION);
    }

    public void setSpecies(ISpecies s) {
        species = s;
    }

    public ISpecies getSpecies() {
        return species;
    }

    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        return (species == null) ? box.getMoleculeList().size() : box.getNMolecules(species);
    }

    public IData getData(IMolecule molecule) {
        if (species != null && molecule.getType() != species) {
            data.x = 0;
            return data;
        }
        double totalKE = 0;
        for (IAtom a : molecule.getChildList()) {
            Vector v = ((IAtomKinetic) a).getVelocity();
            double m = a.getType().getMass();
            totalKE += 0.5 * m * v.squared();
        }
        data.x = (2. / box.getSpace().getD()) * totalKE / molecule.getChildList().size();
        return data;
    }

    public IDataInfo getMoleculeDataInfo() {
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

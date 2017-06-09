/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Quantity;
import etomica.units.Volume;

public class MeterDensitySides extends DataSourceScalar {

    public MeterDensitySides(Box box, ISpecies species, boolean ig) {
        super("density", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1,-1}));
        this.box = box;
        this.species = species;
        this.ig = ig;
    }

    public double getDataAsScalar() {
        IMoleculeList moleculeList = box.getMoleculeList(species);
        int n = 0;
        for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
            Vector p = moleculeList.getMolecule(i).getChildList().getAtom(0).getPosition();
            if (p.getX(0) < 0 == ig) {
                n++;
            }
        }
        return 2*n/box.getBoundary().volume();
    }

    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final ISpecies species;
    protected final boolean ig;
}

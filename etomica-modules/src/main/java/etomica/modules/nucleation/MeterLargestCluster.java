/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.nucleation;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.modules.glass.AtomNbrClusterer;
import etomica.units.dimensions.Quantity;

public class MeterLargestCluster extends DataSourceScalar {

    protected final AtomNbrClusterer clusterer;
    protected final Box box;

    public MeterLargestCluster(Box box) {
        super("largest cluster", Quantity.DIMENSION);
        this.box = box;
        clusterer = new AtomNbrClusterer(box, new AtomTest() {
            @Override
            public boolean test(IAtom a) {
                return true;
            }
        });
        clusterer.setNbrMax(1.5);
    }

    public AtomNbrClusterer getClusterer() {
        return clusterer;
    }

    public double getDataAsScalar() {
        clusterer.findClusters();
        int[] firstAtom = clusterer.getFirstAtom();
        int[] nextAtom = clusterer.getNextAtom();
        int largestCluster = 0;
        for (int i = 0; i < firstAtom.length; i++) {
            if (firstAtom[i] == -1) break;
            int s = 0;
            for (int ii = firstAtom[i]; ii != -1; ii = nextAtom[ii]) {
                s++;
            }
            if (s > largestCluster) largestCluster = s;
        }
        return largestCluster;
    }

}

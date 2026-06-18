/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.nucleation;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramDiscrete;
import etomica.data.types.DataDouble;
import etomica.modules.glass.AtomNbrClusterer;
import etomica.units.dimensions.Quantity;

/**
 * Computes the cluster sizes and returns them as a histogram
 */
public class MeterClusterSizes implements IDataSource {

    protected final AtomNbrClusterer clusterer;
    protected final AccumulatorHistogram histogram;

    public MeterClusterSizes(Box box) {
        histogram = new AccumulatorHistogram(new HistogramDiscrete(1e-9));
        clusterer = new AtomNbrClusterer(box, new AtomTest() {
            @Override
            public boolean test(IAtom a) {
                return true;
            }
        });
        IDataInfo di = new DataDouble.DataInfoDouble("cluster size", Quantity.DIMENSION);
        histogram.putDataInfo(di);
    }

    @Override
    public IData getData() {
        clusterer.findClusters();
        int[] firstAtom = clusterer.getFirstAtom();
        int[] nextAtom = clusterer.getNextAtom();
        histogram.reset();
        DataDouble d = new DataDouble();
        for (int i = 0; i < firstAtom.length; i++) {
            if (firstAtom[i] == -1) break;
            int s = 0;
            for (int ii = firstAtom[i]; ii != -1; ii = nextAtom[ii]) {
                s++;
            }
            d.x = s;
            for (int j = 0; j < s; j++) {
                histogram.putData(d);
            }
        }
        return histogram.getData();
    }

    @Override
    public DataTag getTag() {
        return histogram.getTag();
    }

    @Override
    public IDataInfo getDataInfo() {
        return histogram.getDataInfo();
    }
}

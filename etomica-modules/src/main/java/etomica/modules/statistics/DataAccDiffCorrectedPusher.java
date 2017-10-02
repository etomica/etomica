/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.data.AccumulatorAverage;
import etomica.data.IData;

public class DataAccDiffCorrectedPusher extends DataAccPusher {

    protected final AccumulatorAverage acc;

    public DataAccDiffCorrectedPusher(int idx, DataCollector collector, AccumulatorAverage acc) {
        super(idx, collector);
        this.acc = acc;
    }

    public void putData(IData data) {
        double err = data.getValue(0);
        double cor = data.getValue(1);
        double x = cor <= -0.99999 || cor >= 0.99999 ? Double.NaN : (err * Math.sqrt((1 + cor) / (1 - cor)));
        collector.setData(idx, x * Math.sqrt(acc.getSampleCount()));
    }
}

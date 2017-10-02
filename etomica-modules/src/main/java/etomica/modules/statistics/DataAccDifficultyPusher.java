/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.data.AccumulatorAverage;
import etomica.data.IData;

public class DataAccDifficultyPusher extends DataAccPusher {

    protected final AccumulatorAverage acc;

    public DataAccDifficultyPusher(int idx, DataCollector collector, AccumulatorAverage acc) {
        super(idx, collector);
        this.acc = acc;
    }

    public void putData(IData data) {
        collector.setData(idx, data.getValue(0) * Math.sqrt(acc.getSampleCount()));
    }
}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;

public class DataAccPusher implements IDataSink {

    protected final int idx;
    protected final DataCollector collector;

    public DataAccPusher(int idx, DataCollector collector) {
        this.idx = idx;
        this.collector = collector;
    }

    @Override
    public void putData(IData data) {
        collector.setData(idx, data.getValue(0));
    }

    @Override
    public void putDataInfo(IDataInfo dataInfo) {
        // don't care!
    }

}

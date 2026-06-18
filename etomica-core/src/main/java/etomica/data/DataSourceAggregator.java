/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

public class DataSourceAggregator implements IDataSource {

    protected DataDoubleArray data;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final IDataSource[] dataSources;

    public DataSourceAggregator(IDataSource... dataSources) {
        this.dataSources = dataSources;
        tag = new DataTag();
        setupData();
    }

    protected void setupData() {
        int n = 0;
        Dimension myDimension = null;
        for (IDataSource dataSource : dataSources) {
            n += dataSource.getDataInfo().getLength();
            Dimension d = dataSource.getDataInfo().getDimension();
            if (myDimension == null) myDimension = d;
            if (myDimension != d) myDimension = Null.DIMENSION;
        }
        if (dataInfo != null && dataInfo.getLength() == n) return;
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("aggregated data", myDimension, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);
    }

    @Override
    public IData getData() {
        double[] y = data.getData();
        int idx = 0;
        for (IDataSource dataSource : dataSources) {
            IData data = dataSource.getData();
            int n = data.getLength();
            for (int i = 0; i < n; i++) {
                y[idx + i] = data.getValue(i);
            }
            idx += n;
        }
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        setupData();
        return dataInfo;
    }
}

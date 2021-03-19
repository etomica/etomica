/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataInfoFactory;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;

public class DataProcessorErrorBar extends DataProcessor {
    protected DataFunction data;
    protected String label;

    public DataProcessorErrorBar(String label) {
        this.label = label;
    }

    @Override
    protected IData processData(IData inputData) {
        IData avg = ((DataGroup) inputData).getData(0);
        IData err = ((DataGroup) inputData).getData(1);
        double[] y = data.getData();
        for (int i = 0; i < y.length; i++) {
            y[i] = avg.getValue(i) + err.getValue(i);
        }
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        inputDataInfo = ((DataGroup.DataInfoGroup) inputDataInfo).getSubDataInfo(0);
        data = (DataFunction) inputDataInfo.makeData();
        IDataInfoFactory factory = inputDataInfo.getFactory();
        factory.setLabel(label);
        dataInfo = factory.makeDataInfo();
        dataInfo.addTag(tag);
        return dataInfo;
    }
}

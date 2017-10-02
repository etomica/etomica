/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;

/**
 * DataProcessor that sends on the sum of the incoming data values as a scalar.
 * @author Andrew Schultz
 */
public class DataProcessorSum extends DataProcessor {
    private DataDouble data;
    
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        data = new DataDouble();
        dataInfo = new DataDouble.DataInfoDouble(incomingDataInfo.getLabel(), incomingDataInfo.getDimension());
        return dataInfo;
    }

    public IData processData(IData incomingData) {
        data.x = 0;
        int n = incomingData.getLength();
        for (int i=0; i<n; i++) {
            data.x += incomingData.getValue(i);
        }
        return data;
    }
}

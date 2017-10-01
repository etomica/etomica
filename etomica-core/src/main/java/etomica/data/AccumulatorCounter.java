/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.dimensions.Quantity;

/**
 * Data processor that simply counts the number of times its
 * {@code addData} method is invoked. Output is a DataDouble.
 */
public class AccumulatorCounter extends DataAccumulator {

    public AccumulatorCounter() {
        dataInfo = new DataInfoDouble("Counter", Quantity.DIMENSION);
        data = new DataDouble();
    }
    
    /**
     * Returns null, indicating that any Data type is acceptable for input.
     */
    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        return null;
    }

    /**
     * Does nothing.
     * 
     * @return the DataInfo for the output DataInteger
     */
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        dataInfo.clearTags();
        dataInfo.addTags(incomingDataInfo.getTags());
        dataInfo.addTag(tag);
        return dataInfo;
    }

    /**
     * Increments the counter. Argument is ignored.
     */
    protected boolean addData(IData dummyData) {
        data.x++;
        return true;
    }

    /**
     * Returns the DataInteger with the count.
     */
    public IData getData() {
        return data;
    }

    /**
     * Sets count to zero.
     */
    public void reset() {
        data.x = 0;
    }

    private static final long serialVersionUID = 1L;

    protected final DataDouble data;
}

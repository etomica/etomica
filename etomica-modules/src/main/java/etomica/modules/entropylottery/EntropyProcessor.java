/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Calculates the entropy of an incoming probability density
 * using Stirling's approximation.
 * @author Andrew Schultz
 */
public class EntropyProcessor extends DataProcessor {

    public EntropyProcessor() {
        super();
    }
    
    public DataPipe getDataCaster(IEtomicaDataInfo newDataInfo) {
        // we actually just want DataArithmetic, but array is probably
        // reasonable
        if (!(newDataInfo instanceof DataInfoDoubleArray)) {
            return new CastToDoubleArray();
        }
        return null;
    }
    
    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("entropy", Null.DIMENSION);
        dataInfo.addTags(incomingDataInfo.getTags());
        dataInfo.addTag(tag);
        return dataInfo;
    }
        

    public IData processData(IData incomingData) {
        double sum = 0;
        DataDoubleArray inData = (DataDoubleArray)incomingData;
        double totalCount = 0;
        for (int i=0; i<incomingData.getLength(); i++) {
            totalCount += inData.getValue(i);
        }
        // - k sum (n log n/N)
        for (int i=0; i<incomingData.getLength(); i++) {
            double x = inData.getValue(i);
            if (x > 0) {
                sum += x * Math.log(x/totalCount);
            }
        }
        data.x = -sum;
        return data;
    }

    private static final long serialVersionUID = 1L;
    protected DataDouble data;
}

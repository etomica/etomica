/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.dimensions.Quantity;

/**
 * Takes output from MeterChainLength as input and outputs the number-average
 * molecular weight.
 * 
 * @author Andrew Schultz
 */
public class MolecularWeightAvg extends DataProcessor {

    public MolecularWeightAvg() {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Avg MW", Quantity.DIMENSION);
    }

    protected IData processData(IData inputData) {
        double sum = 0, sum2 = 0;
        for (int i=0; i<inputData.getLength(); i++) {
            double v = inputData.getValue(i);
            sum += v;
            sum2 += v / (i+1);
        }
        data.x = sum / sum2;
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return dataInfo;
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
}

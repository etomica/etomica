/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.surfacetension;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.space.Tensor;
import etomica.units.dimensions.Null;

public class DataProcessorSurfaceTension extends DataProcessor {

    protected final DataInfoDouble dataInfo = new DataInfoDouble("surface tension", Null.DIMENSION);
    protected final DataDouble data = new DataDouble();

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        if (!(inputDataInfo instanceof DataInfoTensor)) throw new RuntimeException("I want a DataTensor");
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        Tensor t = ((DataTensor)inputData).x;
        data.x = 0.5*(t.component(0,0) - 0.5*(t.component(1,1) + t.component(2,2)));
        return data;
    }
}

package etomica.surfacetension;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataTensor;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.space.Tensor;
import etomica.units.dimensions.Null;

public class DataProcessorSurfaceTension extends DataProcessor {
    DataInfoDouble dataInfo = new DataInfoDouble("surface tension", Null.DIMENSION);

    DataDouble data = new DataDouble();

    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        if (!(inputDataInfo instanceof DataInfoTensor)) throw new RuntimeException("I want a DataTensor");
        return null;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        Tensor t = ((DataTensor)inputData).x;
        data.x = 0.5*(t.component(0,0) - 0.5*(t.component(1,1) + t.component(2,2)));
        return data;
    }
}

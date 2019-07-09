package etomica.cavity;

import etomica.box.Box;
import etomica.data.DataProcessorForked;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;

/**
 * DataProcessor that determines g(sigma) from P (measured with the virial)
 */
public class DataProcessorPContactG extends DataProcessorForked {

    protected final DataGroup data = new DataGroup(new DataDouble[]{new DataDouble(), new DataDouble(), new DataDouble()});
    protected final Box box;

    public DataProcessorPContactG(Box box) {
        this.box = box;
        DataDouble.DataInfoDouble subDataInfo = new DataDouble.DataInfoDouble("stuff", Null.DIMENSION);
        dataInfo = new DataGroup.DataInfoGroup("pressure", Pressure.DIMENSION, new IDataInfo[]{subDataInfo, subDataInfo, subDataInfo});
        dataInfo.addTag(tag);
    }

    @Override
    protected IData processData(IData inputData) {
        double rho = box.getLeafList().size() / box.getBoundary().volume();
        double gsigma = inputData.getValue(0);
        // g(sigma) = (p/rho - 1) 1.5 / (pi rho)
        // p = rho + (2/3) pi rho^2 g(sigma)
        ((DataDouble) data.getData(0)).x = rho * (1 + 2.0 / 3.0 * Math.PI * rho * gsigma);
        gsigma = inputData.getValue(1);
        ((DataDouble) data.getData(1)).x = rho * (1 + 2.0 / 3.0 * Math.PI * rho * gsigma);
        double e = inputData.getValue(2);
        ((DataDouble) data.getData(2)).x = rho * 2.0 / 3.0 * Math.PI * rho * e;
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return dataInfo;
    }
}

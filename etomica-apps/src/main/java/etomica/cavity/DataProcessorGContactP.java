package etomica.cavity;

import etomica.box.Box;
import etomica.data.DataProcessorForked;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Null;

/**
 * DataProcessor that determines g(sigma) from P (measured with the virial)
 */
public class DataProcessorGContactP extends DataProcessorForked {

    protected final DataGroup data = new DataGroup(new DataDouble[]{new DataDouble(), new DataDouble(), new DataDouble()});
    protected final Box box;

    public DataProcessorGContactP(Box box) {
        this.box = box;
        DataDouble.DataInfoDouble subDataInfo = new DataDouble.DataInfoDouble("stuff", Null.DIMENSION);
        dataInfo = new DataGroup.DataInfoGroup("g(sigma)", Energy.DIMENSION, new IDataInfo[]{subDataInfo, subDataInfo, subDataInfo});
        dataInfo.addTag(tag);
    }

    @Override
    protected IData processData(IData inputData) {
        double rho = box.getLeafList().size() / box.getBoundary().volume();
        double p = inputData.getValue(0);
        ((DataDouble) data.getData(0)).x = (p / rho - 1) * 1.5 / Math.PI / rho;
        p = inputData.getValue(1);
        ((DataDouble) data.getData(1)).x = (p / rho - 1) * 1.5 / Math.PI / rho;
        double e = inputData.getValue(2);
        ((DataDouble) data.getData(2)).x = (e / rho) * 1.5 / Math.PI / rho;
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return dataInfo;
    }
}

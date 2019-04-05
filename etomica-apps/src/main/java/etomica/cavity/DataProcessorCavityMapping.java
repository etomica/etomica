package etomica.cavity;

import etomica.data.DataProcessor;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataFunction;

public class DataProcessorCavityMapping extends DataProcessor {

    protected final P2HardSphereCavity p2;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag tag;

    public DataProcessorCavityMapping(P2HardSphereCavity p2) {
        super();
        this.p2 = p2;
        tag = new DataTag();
    }

    @Override
    protected IData processData(IData inputData) {
        int l = inputData.getLength();
        double[] y = data.getData();
        double dx = p2.getCollisionDiameter() / l;
        for (int i = 0; i < l; i++) {
            double x = inputData.getValue(i);
            for (int j = 0; j < i; j++) {
                y[j] += x;
            }
        }
        for (int i = 0; i < l; i++) {
            y[i] /= 4 * Math.PI;
        }
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = (DataFunction.DataInfoFunction) inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = (DataFunction) inputDataInfo.makeData();
        return dataInfo;
    }
}

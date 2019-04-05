package etomica.cavity;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorHard;

/**
 * This processor scales up the measured cavity function so that its contact
 * value is equal to the RDF contact value.
 */
public class DataProcessorCavity extends DataProcessor {

    protected DataFunction data;
    protected final IntegratorHard integrator;
    protected final P2HardSphereCavity p2;

    public DataProcessorCavity(IntegratorHard integrator, P2HardSphereCavity p2) {
        super();
        this.integrator = integrator;
        this.p2 = p2;
    }

    @Override
    protected IData processData(IData inputData) {
        long totalCollision = integrator.getCollisionCount();
        long internalCollision = p2.getInternalCount();
        double fac = 1;
        if (internalCollision > 0) {
            long externalCollision = totalCollision - internalCollision;
            fac = externalCollision / (double) internalCollision;
        }
        DataDoubleArray rData = ((DataFunction.DataInfoFunction) dataInfo).getXDataSource().getIndependentData(0);
        double[] y = data.getData();
        for (int i = 0; i < inputData.getLength(); i++) {
            if (rData.getValue(i) > p2.getCollisionDiameter()) y[i] = 0;
            else y[i] = inputData.getValue(i) * fac;
        }
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = (DataFunction) inputDataInfo.makeData();
        return dataInfo;
    }
}

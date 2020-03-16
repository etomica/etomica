package etomica.cavity;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataInfoFactory;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorHard;

/**
 * This processor scales up the measured cavity function so that its contact
 * value is equal to the RDF contact value.
 */
public class DataProcessorCavity extends DataProcessor implements IntegratorHard.CollisionListener {

    protected DataFunction data;
    protected long internalCollisions, totalCollisions;
    protected double sigma;

    public DataProcessorCavity(IntegratorHard integrator) {
        super();
        integrator.addCollisionListener(this);
    }

    public void reset() {
        totalCollisions = internalCollisions = 0;
    }

    public void collisionAction(IntegratorHard.Agent agent) {
        totalCollisions++;
        P2HardSphereCavity p2 = (P2HardSphereCavity) agent.collisionPotential;
        P2HardSphereCavity.CollisionType cType = p2.getLastCollisionType();
        sigma = p2.getCollisionDiameter();
        if (cType == P2HardSphereCavity.CollisionType.ESCAPE) {
            internalCollisions++;
        } else if (cType == P2HardSphereCavity.CollisionType.INTERNAL_BOUNCE) {
            internalCollisions++;
        }
    }

    @Override
    protected IData processData(IData inputData) {
        double fac = 1;
        if (internalCollisions > 0) {
            long externalCollision = totalCollisions - internalCollisions;
            fac = externalCollision / (double) internalCollisions;
        }
        DataDoubleArray rData = ((DataFunction.DataInfoFunction) dataInfo).getXDataSource().getIndependentData(0);
        double[] y = data.getData();
        for (int i = 0; i < inputData.getLength(); i++) {
            if (rData.getValue(i) > sigma) y[i] = 0;
            else y[i] = inputData.getValue(i) * fac;
        }
        reset();
        return data;
    }

    @Override
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        IDataInfoFactory factory = inputDataInfo.getFactory();
        factory.setLabel("y(r)");
        dataInfo = factory.makeDataInfo();
        dataInfo.addTag(tag);
        data = (DataFunction) inputDataInfo.makeData();
        return dataInfo;
    }
}

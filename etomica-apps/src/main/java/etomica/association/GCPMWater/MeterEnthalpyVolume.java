package etomica.association.GCPMWater;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Energy;

public class MeterEnthalpyVolume implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final IntegratorBox IB;
    protected final double pressure;

    public MeterEnthalpyVolume(IntegratorBox IB, double pressure) {
        this.IB = IB;
        this.pressure = pressure;
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("energies!", Energy.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    @Override
    public IData getData() {
        double PE = IB.getPotentialEnergy();
        double Volume = IB.getBox().getBoundary().volume();
        double[] x = data.getData();
        x[0] = PE + (pressure*Volume);
        x[1] = Volume;
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
}

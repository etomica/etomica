package etomica.modules.interfacial;

import etomica.api.IData;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Quantity;

public class MeterDensityProfileFromForceSum implements IEtomicaDataSource, DataSourceIndependent {

    public MeterDensityProfileFromForceSum(PotentialCalculationForcePressureBinSum forceSum) {
        this.forceSum = forceSum;
        data = new DataFunction(new int[]{0});
        setupData();
        tag = new DataTag();
    }
    
    protected void setupData() {
        double[] densityProfile = forceSum.getDensityProfile();
        if (densityProfile != data.getData()) {
            int D = forceSum.space.D();
            xData = new DataDoubleArray(densityProfile.length);
            double binSize = forceSum.binSize;
            double x = -forceSum.halfL + 0.5*binSize;
            double[] xValues = xData.getData();
            for (int i=0; i<densityProfile.length; i++) {
                xValues[i] = x;
                x += binSize;
            }
            xDataInfo = new DataInfoDoubleArray("x position", Length.DIMENSION, new int[]{densityProfile.length});
            data = new DataFunction(new int[]{densityProfile.length}, densityProfile);
            dataInfo = new DataInfoFunction("density profile", new CompoundDimension(
                new Dimension[]{Quantity.DIMENSION, Length.DIMENSION}, new double[]{1, -D}), this);
        }
    }
    
    public IData getData() {
        if (data.getData() != forceSum.getDensityProfile()) {
            setupData();
        }
        return data;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final PotentialCalculationForcePressureBinSum forceSum;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag tag;
}

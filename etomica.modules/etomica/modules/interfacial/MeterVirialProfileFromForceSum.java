package etomica.modules.interfacial;

import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.space.Tensor;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Pressure;
import etomica.units.Pressure2D;

public class MeterVirialProfileFromForceSum implements DataSource, DataSourceIndependent {

    public MeterVirialProfileFromForceSum(PotentialCalculationForcePressureBinSum forceSum) {
        this.forceSum = forceSum;
        data = new DataGroup(new Data[]{new DataFunction(new int[]{0})});
        setupData();
        tag = new DataTag();
    }
    
    protected void setupData() {
        Tensor[] virialProfile = forceSum.getPressureTensorProfile();

        if (virialProfile.length != ((DataFunction)data.getData(0)).getLength()) {
            int D = forceSum.space.D();
            int len = virialProfile.length;
            xData = new DataDoubleArray(len);
            double binSize = forceSum.binSize;
            double x = -forceSum.halfL + 0.5*binSize;
            double[] xValues = xData.getData();
            for (int i=0; i<len; i++) {
                xValues[i] = x;
                x += binSize;
            }
            xDataInfo = new DataInfoDoubleArray("x position", Length.DIMENSION, new int[]{len});

            DataFunction[] dataFunctions = new DataFunction[D];
            Dimension pDim = Pressure.DIMENSION;
            if (D == 2) {
                pDim = Pressure2D.DIMENSION;
            }
            DataInfoFunction dataInfoFunction = new DataInfoFunction("virial profile", pDim, this);
            DataInfoFunction[] dataInfoFunctions = new DataInfoFunction[D];
            y = new double[D][0];
            for (int i=0; i<D; i++) {
                dataFunctions[i] = new DataFunction(new int[]{len});
                y[i] = dataFunctions[i].getData();
                dataInfoFunctions[i] = dataInfoFunction;
            }
            data = new DataGroup(dataFunctions);
            dataInfo = new DataInfoGroup("virial profiles", pDim, dataInfoFunctions);

        }
    }
    
    public Data getData() {
        setupData();

        Tensor[] virialProfile = forceSum.getPressureTensorProfile();
        int D = virialProfile[0].D();
        for (int i=0; i<virialProfile.length; i++) {
            for (int j=0; j<D; j++) {
                y[j][i] = virialProfile[i].component(j,j);
            }
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

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final PotentialCalculationForcePressureBinSum forceSum;
    protected DataGroup data;
    protected DataInfoGroup dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag tag;
    protected double[][] y;
}

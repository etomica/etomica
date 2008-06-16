package etomica.modules.interfacial;

import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.IDataInfoFactory;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.space.ISpace;
import etomica.units.Area;
import etomica.units.DimensionRatio;
import etomica.units.Energy;
import etomica.units.Length;

public class DataProcessorInterfacialTensionProfile extends DataProcessor {

    public DataProcessorInterfacialTensionProfile(PotentialCalculationForcePressureBinSum forceSum) {
        this.forceSum = forceSum;
        virialData = new double[forceSum.space.D()][0];
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }

    public IBox getBox() {
        return box;
    }

    protected Data processData(Data inputData) {
        DataGroup dataGroup = (DataGroup)inputData;
        int D = virialData.length;
        for (int i=0; i<D; i++) {
            virialData[i] = ((DataFunction)dataGroup.getData(i)).getData();
        }
        int nBins = data.getArrayShape(0);
        double[] tension = data.getData();
        for (int i=0; i<nBins; i++) {
            tension[i] += (D-1)*virialData[0][i];
        }
        for (int j=1; j<D; j++) {
            for (int i=0; i<nBins; i++) {
                tension[i] -= virialData[j][i];
            }
        }

        double area = 1;
        IVector dim = box.getBoundary().getDimensions();
        for (int i=1; i<dim.getD(); i++) {
            area *= dim.x(i);
        }
        double fac = 0.25/area/forceSum.binSize;
        for (int i=0; i<nBins; i++) {
            tension[i] *= fac;
        }
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        DataInfoFunction dataInfo0 = (DataInfoFunction)((DataInfoGroup)inputDataInfo).getSubDataInfo(0);
        data = (DataFunction)dataInfo0.makeData();
        IDataInfoFactory dataInfoFactory = dataInfo0.getFactory();
        dataInfoFactory.setDimension(new DimensionRatio(Energy.DIMENSION, ((DataInfoGroup)inputDataInfo).getNDataInfo() == 2 ? Length.DIMENSION : Area.DIMENSION));
        dataInfoFactory.setLabel("Interfacial tension profile");
        return dataInfoFactory.makeDataInfo();
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    protected DataFunction data;
    protected final double[][] virialData;
    protected IBox box;
    protected final PotentialCalculationForcePressureBinSum forceSum;
}

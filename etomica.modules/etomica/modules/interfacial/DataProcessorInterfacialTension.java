package etomica.modules.interfacial;

import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataTensor;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Area;
import etomica.units.DimensionRatio;
import etomica.units.Energy;
import etomica.units.Length;

public class DataProcessorInterfacialTension extends DataProcessor {

    public DataProcessorInterfacialTension(ISpace space) {
        this.space = space;
        data = new DataDouble();
    }

    public void setBox(IBox newBox) {
        box = newBox;
    }
    
    public IBox getBox() {
        return box;
    }

    public Data processData() {
        return data;
    }

    protected Data processData(Data inputData) {
        double area = 1;
        IVector dim = box.getBoundary().getDimensions();
        for (int i=1; i<dim.getD(); i++) {
            area *= dim.x(i);
        }
        data.x = 0;
        double vxx = inputData.getValue(0);
        for (int i=1; i<inputData.getLength(); i++) {
            data.x += vxx - inputData.getValue(i);
        }
        data.x /= 4*area;
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return new DataInfoDouble("Interfacial tension", new DimensionRatio(Energy.DIMENSION, space.D() == 2 ? Length.DIMENSION : Area.DIMENSION));
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    protected final ISpace space;
    protected IBox box;
    protected final DataDouble data;
}

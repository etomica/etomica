package etomica.modules.interfacial;

import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.space.ISpace;
import etomica.units.Area;
import etomica.units.DimensionRatio;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Data Processor that takes the virial components as input data and returns
 * the interfacial tension.  This class must be told which dimension the
 * surfaces exist in.
 *
 * @author Andrew Schultz
 */
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

    /**
     * Sets the dimension of the surfaces.  0 (x) means that the surfaces exist
     * in the yz plane.
     */
    public void setSurfaceDimension(int newDim) {
        surfaceDim = newDim;
    }

    /**
     * Returns the dimension of the surfaces.  0 (x) means that the surfaces
     * exist in the yz plane.
     */
    public int getSurfaceDimension() {
        return surfaceDim;
    }

    protected Data processData(Data inputData) {
        double area = 1;
        IVector dim = box.getBoundary().getDimensions();
        int D = dim.getD();
        for (int i=0; i<D; i++) {
            if (i == surfaceDim) continue;
            area *= dim.x(i);
        }
        data.x = 0;
        data.x = inputData.getValue(0);
        for (int i=0; i<D; i++) {
            if (i == surfaceDim) continue;
            data.x -= inputData.getValue(i)/(D-1);
        }
        data.x /= 2*area;
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
    protected int surfaceDim;
}

package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.IData;
import etomica.api.IMoleculeList;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.Energy;

/**
 * Meter that calculates the Boltzmann-factored harmonic energy of each normal mode for a 
 * configuration given eigenvectors and omegas corresponding to wave vectors.
 * 
 * @author Andrew Schultz
 */
public class MeterHarmonicCoordinate implements IEtomicaDataSource {

    public MeterHarmonicCoordinate(CoordinateDefinition coordinateDefinition, NormalModes normalModes, int[] modes) {
        this.coordinateDefinition = coordinateDefinition;
        this.normalModes = normalModes;
        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{0});
        tag = new DataTag();
        this.modes = modes;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    

    public IData getData() {
        BasisCell cell = coordinateDefinition.getBasisCells()[0];
        int coordinateDim = coordinateDefinition.getCoordinateDim();

        IMoleculeList molecules = cell.molecules;
        double[] u = coordinateDefinition.calcU(molecules);
        double sqrtCells = Math.sqrt(coordinateDefinition.getBasisCells().length);

        double[] x = data.getData();
        
        for (int i=0; i<modes.length; i++) {
            double realCoord = 0;
            for (int j=0; j<coordinateDim; j++) {
                realCoord += u[j] * eigenvectors[0][modes[i]][j] * sqrtCells;
            }
            x[i] = realCoord;
        }
        return data;
    }
    
    public IBox getBox() {
        return coordinateDefinition.getBox();
    }

    public void setBox(IBox newBox) {
        normalModes.getWaveVectorFactory().makeWaveVectors(newBox);
        setEigenvectors(normalModes.getEigenvectors(newBox));

        dataInfo = new DataInfoDoubleArray("Harmonic Coordinates", Energy.DIMENSION, new int[]{modes.length});
        data = new DataDoubleArray(new int[]{modes.length});
    }
    
    public void setEigenvectors(double[][][] eigenVectors) {
        this.eigenvectors = (double[][][])eigenVectors.clone();
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    private static final long serialVersionUID = 1L;
    protected final CoordinateDefinition coordinateDefinition;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    private final DataTag tag;
    protected double[][][] eigenvectors;
    protected String name;
    protected NormalModes normalModes;
    protected int[] modes;
}

package etomica.normalmode;

import etomica.atom.AtomSet;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.phase.Phase;
import etomica.units.Energy;

/**
 * Meter that calculates the Boltzmann-factored harmonic energy of each normal mode for a 
 * configuration given eigenvectors and omegas corresponding to wave vectors.
 * 
 * @author Andrew Schultz
 */
public class MeterHarmonicCoordinate implements DataSource {

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

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    

    public Data getData() {
        BasisCell cell = coordinateDefinition.getBasisCells()[0];
        int coordinateDim = coordinateDefinition.getCoordinateDim();

        AtomSet molecules = cell.molecules;
        double[] u = coordinateDefinition.calcU(molecules);
        double sqrtCells = Math.sqrt(coordinateDefinition.getBasisCells().length);
        // divide by sqrt2.  Don't ask why.
        double sqrt2 = Math.sqrt(2);
        for (int i = 0; i < coordinateDim; i++) {
            u[i] *= sqrtCells / sqrt2;
        }

        double[] x = data.getData();
        
        for (int i=0; i<modes.length; i++) {
            double realCoord = 0;
            for (int j=0; j<coordinateDim; j++) {
                realCoord += u[j] * eigenvectors[0][j][i];
            }
            x[i] = realCoord;
        }
        return data;
    }
    
    public Phase getPhase() {
        return coordinateDefinition.getPhase();
    }

    public void setPhase(Phase newPhase) {
        normalModes.getWaveVectorFactory().makeWaveVectors(newPhase);
        setEigenvectors(normalModes.getEigenvectors(newPhase));

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

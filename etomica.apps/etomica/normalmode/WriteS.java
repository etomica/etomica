package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.Action;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.IVector;

/**
 * Class that writes out S from the MeterNormalMode, calculates
 * eigenvectors/values (writes them out) and calculates the harmonic free
 * energy from them.
 *
 * @author Andrew Schultz
 */
public class WriteS implements Action {

    public WriteS() {
        temperature = 1.0;
        doOverwrite = true;
    }
    
    public void setMeter(MeterNormalMode meter) {
        meterNormalMode = meter;
    }
    
    public void setFilename(String newFilename) {
        filename = newFilename;
        count = 0;
    }
    
    public void setWaveVectorFactory(WaveVectorFactory newWaveVectorFactory) {
        waveVectorFactory = newWaveVectorFactory;
    }
    
    public void setOverwrite(boolean newDoOverwrite) {
        doOverwrite = newDoOverwrite;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void actionPerformed() {
        // normalize averages
        DataGroup normalModeData = (DataGroup) meterNormalMode.getData();
        int callCount = meterNormalMode.getCallCount();

        // write wave vectors (to filename.k) and simulation results (to
        // filename.S) to file
        IVector[] waveVectors = waveVectorFactory.getWaveVectors();
        double[] coefficients = waveVectorFactory.getCoefficients();

        String thisFilename = filename;
        
        try {
            int coordinateDim = meterNormalMode.getCoordinateDefinition()
                    .getCoordinateDim();
            if (!doOverwrite) {
                thisFilename += "_" + count;
                count++;
            }
            FileWriter fileWriterK = new FileWriter(thisFilename + ".k");
            FileWriter fileWriterS = new FileWriter(thisFilename + ".S");
            for (int k = 0; k < waveVectors.length; k++) {
                // write the wavevector with its coefficient
                fileWriterK.write(Double.toString(coefficients[k]));
                for (int j = 0; j < waveVectors[k].getD(); j++) {
                    fileWriterK.write(" " + waveVectors[k].x(j));
                }
                fileWriterK.write("\n");

                // write the (coordDim x coordDim) S array for the current
                // wavevector
                DataDoubleArray dataS = (DataDoubleArray) normalModeData.getData(k);
                for (int j = 0; j < coordinateDim; j++) {
                    fileWriterS.write(Double.toString(dataS.getValue(j * coordinateDim)/callCount));
                    for (int l = 1; l < coordinateDim; l++) {
                        fileWriterS.write(" " + dataS.getValue(j * coordinateDim + l)/callCount);
                    }
                    fileWriterS.write("\n");
                }
            }
            fileWriterK.close();
            fileWriterS.close();
        } catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data " + e);
        }
        
        NormalModeEigenGetter.doit(thisFilename);

        BasisCell[] cells = meterNormalMode.getCoordinateDefinition().getBasisCells();
        CalcHarmonicA.doit(thisFilename, meterNormalMode.getPhase().getSpace().D(), 1.0, temperature, cells[0].molecules.getAtomCount(), cells.length);
    }

    protected MeterNormalMode meterNormalMode;
    protected int count;
    protected boolean doOverwrite;
    protected String filename;
    protected WaveVectorFactory waveVectorFactory;
    protected double temperature;
}

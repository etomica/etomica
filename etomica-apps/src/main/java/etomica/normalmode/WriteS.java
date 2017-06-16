/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.IAction;
import etomica.space.Vector;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.space.Space;

/**
 * Class that writes out S from the MeterNormalMode, calculates
 * eigenvectors/values (writes them out) and calculates the harmonic free
 * energy from them.
 *
 * @author Andrew Schultz
 */
public class WriteS implements IAction {

    public WriteS(Space _space) {
        temperature = 1.0;
        doOverwrite = true;
        space = _space;
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
        Vector[] waveVectors = waveVectorFactory.getWaveVectors();
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
                    fileWriterK.write(" " + waveVectors[k].getX(j));
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
        
        NormalModeEigenGetter.doit(thisFilename, space.D());

        int numMolecules = meterNormalMode.getBox().getMoleculeList().getMoleculeCount();
        NormalModesFromFile normalModes = new NormalModesFromFile(thisFilename, space.D());
        normalModes.setTemperature(temperature);
        lastA = CalcHarmonicA.doit(normalModes, space.D(), temperature, numMolecules);
    }
    
    public double getLastA() {
        return lastA;
    }

    protected MeterNormalMode meterNormalMode;
    protected int count;
    protected boolean doOverwrite;
    protected String filename;
    protected WaveVectorFactory waveVectorFactory;
    protected double temperature;
    private final Space space;
    protected double lastA;
}

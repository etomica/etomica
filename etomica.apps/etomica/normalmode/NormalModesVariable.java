package etomica.normalmode;

import etomica.api.IVectorMutable;
import etomica.space.ISpace;

public class NormalModesVariable {

    public NormalModesVariable(ISpace space, int nModes, CoordinateDefinition coordinateDefinition) {
        eigenVectors = new double[nModes][coordinateDefinition.getCoordinateDim()];
        waveVectors = new IVectorMutable[nModes];
        for (int i=0; i<nModes; i++) {
            waveVectors[i] = space.makeVector();
        }
        phaseAngles = new double[nModes];
    }
    
    public double[][] getEigenVectors() {
        return eigenVectors;
    }

    public IVectorMutable[] getWaveVectors() {
        return waveVectors;
    }
    
    public double[] getPhaseAngles() {
        return phaseAngles;
    }

    protected IVectorMutable[] waveVectors;
    protected double[][] eigenVectors;
    protected double[] phaseAngles;
}

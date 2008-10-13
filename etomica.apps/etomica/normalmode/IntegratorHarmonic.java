package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.integrator.IntegratorMD;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.ISpace;

/**
 * 
 * Integrator that is used for the harmonic osciallating system
 * 
 * 
 * @author Tai Boon Tan
 */
public class IntegratorHarmonic extends IntegratorMD {

    public IntegratorHarmonic(IRandom random, double timeStep, double temperature, ISpace _space) {
        super(null,random, timeStep, temperature, _space);
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setOmegaSquared(double[][] omega2, double[] coeff) {
        omega = new double[omega2.length][omega2[0].length];
        for (int i=0; i<omega.length; i++) {
            for (int j=0; j<omega[i].length; j++) {
                omega[i][j] = Math.sqrt(omega2[i][j]);
            }
        }
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void setWaveVectors(IVector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }
    
    public void setWaveVectorCoefficients(double[] newWaveVectorCoefficients) {
        waveVectorCoefficients = newWaveVectorCoefficients;
    }
    
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setBox(IBox newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        u = new double[coordinateDim];

        rRand = new double[waveVectors.length][coordinateDim];
        iRand = new double[waveVectors.length][coordinateDim];
    }

    public void doStepInternal() {
    	super.doStepInternal();
        iterator.reset();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();

        double sqrtT = Math.sqrt(temperature);

        for (int iVector=0; iVector<waveVectors.length; iVector++) {
            for (int j=0; j<coordinateDim; j++) {
                if (omega[iVector][j] == Double.POSITIVE_INFINITY) continue;
                
                //generate real and imaginary parts of random normal-mode coordinate Q
                rRand[iVector][j] = Math.cos(omega[iVector][j]*currentTime) *sqrtT;
                iRand[iVector][j] = Math.sin(omega[iVector][j]*currentTime) *sqrtT;
                
            }
        }
        
        
        for (int iCell = 0; iCell<cells.length; iCell++) {

            BasisCell cell = cells[iCell];
            for (int i=0; i<coordinateDim; i++) {
                u[i] = 0;
            }
            //loop over wavevectors and sum contribution of each to the generalized coordinates
            for (int iVector=0; iVector<waveVectors.length; iVector++) {
                double kR = waveVectors[iVector].dot(cell.cellPosition);//getLatticePositions()[atomCount]);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                
                for (int i=0; i<coordinateDim; i++) {
                    for (int j=0; j<coordinateDim; j++) {
                        u[j] += waveVectorCoefficients[iVector]*eigenVectors[iVector][i][j]*
                                  2.0*(rRand[iVector][i]*coskR - iRand[iVector][i]*sinkR);
                    }
                }
            }
            double normalization = 1/Math.sqrt(cells.length);
            for (int i=0; i<coordinateDim; i++) {
                u[i] *= normalization;
            }
            coordinateDefinition.setToU(cell.molecules, u);
        }
    }
    



    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorAllMolecules iterator;
    private double[][] omega;
    private double[][][] eigenVectors;
    private IVector[] waveVectors;
    private double[] waveVectorCoefficients;
    protected double[] u;
    protected double[][] rRand;
    protected double[][] iRand;
    protected double lastEnergy;
    protected double temperature;
    protected boolean isRejectable;
    protected double[][] uOld;
}

package etomica.normalmode;

import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.integrator.mcmove.MCMoveTracker;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.util.IRandom;

public class MCMoveHarmonic extends MCMovePhase {

    public MCMoveHarmonic(Simulation sim) {
        this(sim.getPotentialMaster(), sim.getRandom());
    }
    
    public MCMoveHarmonic(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster, new MCMoveTracker());
        this.random = random;
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public void setOmegaSquared(double[][] omega2, double[] coeff) {
        stdDev = new double[omega2.length][omega2[0].length];
        for (int i=0; i<stdDev.length; i++) {
            for (int j=0; j<stdDev[i].length; j++) {
                stdDev[i][j] = Math.sqrt(1.0/(2.0*omega2[i][j]*coeff[i]));
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
    
    public void setPhase(Phase newPhase) {
        super.setPhase(newPhase);
        iterator.setPhase(newPhase);

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        u = new double[coordinateDim];

        rRand = new double[waveVectors.length][coordinateDim];
        iRand = new double[waveVectors.length][coordinateDim];
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public boolean doTrial() {
        iterator.reset();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();

        lastEnergy = 0;
        double sqrtT = Math.sqrt(temperature);

        for (int iVector=0; iVector<waveVectors.length; iVector++) {
            for (int j=0; j<coordinateDim; j++) {
                if (stdDev[iVector][j] == 0) continue;
                //generate real and imaginary parts of random normal-mode coordinate Q
                double realGauss = random.nextGaussian();
                double imaginaryGauss = random.nextGaussian();
                rRand[iVector][j] = realGauss * stdDev[iVector][j] * sqrtT;
                iRand[iVector][j] = imaginaryGauss * stdDev[iVector][j] * sqrtT;
                //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
                // it's hard to know which.  So long as we don't put an atom at the origin
                // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
                if (waveVectorCoefficients[iVector] == 0.5) imaginaryGauss = 0;
                lastEnergy += 0.5 * (realGauss*realGauss + imaginaryGauss*imaginaryGauss);
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
        return true;
    }
    
    public double getA() {
        // return 1 to guarantee success
        return 1;
    }

    public double getB() {
        // return 0 to guarantee success
        return 0;
    }
    
    /**
     * Returns the harmonic energy of the configuration based on the last
     * harmonic move made by this MC Move.
     */
    public double getLastTotalEnergy() {
        return lastEnergy;
    }
    
    public void acceptNotify() {
    }

    public double energyChange() {
        return 0;
    }

    public void rejectNotify() {
        throw new RuntimeException("This move should never be rejected");
    }

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorAllMolecules iterator;
    private double[][] stdDev;
    private double[][][] eigenVectors;
    private IVector[] waveVectors;
    private double[] waveVectorCoefficients;
    protected double[] u;
    protected double[][] rRand;
    protected double[][] iRand;
    protected final IRandom random;
    protected double lastEnergy;
    protected double temperature;
}

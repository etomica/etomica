package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.integrator.mcmove.MCMoveTracker;
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

        coordinateDefinition.setPhase(newPhase);
//        latticePositions = coordinateDefinition.getLatticePositions();
//
//        iterator.reset();
//        int atomCount = 0;
//        while (iterator.hasNext()) {
//            latticePositions[atomCount] = phase.getSpace().makeVector();
//            Atom atom = iterator.nextAtom();
//            latticePositions[atomCount].E(atom.getType().getPositionDefinition().position(atom));
//            atomCount++;
//        }

        int coordinateDim = coordinateDefinition.getCoordinateDim();
        u = new double[coordinateDim];

        rRand = new double[waveVectors.length][coordinateDim];
        iRand = new double[waveVectors.length][coordinateDim];
        
        normalization = 1/Math.sqrt(phase.getSpeciesMaster().moleculeCount());
        
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public boolean doTrial() {
        iterator.reset();
        int coordinateDim = coordinateDefinition.getCoordinateDim();

        for (int iVector=0; iVector<waveVectors.length; iVector++) {
            for (int j=0; j<coordinateDim; j++) {
                //generate real and imaginary parts of random normal-mode coordinate Q
                rRand[iVector][j] = random.nextGaussian() * stdDev[iVector][j];
                iRand[iVector][j] = random.nextGaussian() * stdDev[iVector][j];
            }
        }
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            for (int i=0; i<coordinateDim; i++) {
                u[i] = 0;
            }
            //loop over wavevectors and sum contribution of each to the generalized coordinates
            for (int iVector=0; iVector<waveVectors.length; iVector++) {
                double kR = waveVectors[iVector].dot(coordinateDefinition.getLatticePosition(atom));//getLatticePositions()[atomCount]);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                
                for (int i=0; i<coordinateDim; i++) {
                    for (int j=0; j<coordinateDim; j++) {
                        u[j] += waveVectorCoefficients[iVector]*eigenVectors[iVector][i][j]*
                                  2.0*(rRand[iVector][i]*coskR - iRand[iVector][i]*sinkR);
                    }
                }
            }
            for (int i=0; i<coordinateDim; i++) {
                u[i] *= normalization;
            }
            coordinateDefinition.setToU(atom, u);
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
    protected double normalization;
    protected final IRandom random;
}

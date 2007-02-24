package etomica.normalmode;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.integrator.mcmove.MCMoveTracker;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;

public class MCMoveHarmonic extends MCMovePhase {

    public MCMoveHarmonic(PotentialMaster potentialMaster) {
        super(potentialMaster, new MCMoveTracker());
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setNormalCoordWrapper(NormalCoordMapper newNormalCoordWrapper) {
        normalCoordMapper = newNormalCoordWrapper;
        normalDim = normalCoordMapper.getNormalDim();
    }
    
    public NormalCoordMapper getNormalCoordWrapper() {
        return normalCoordMapper;
    }

    public void setEigenValues(double[][] newEigenValues) {
        eigenValuesSqrt = new double[newEigenValues.length][newEigenValues[0].length];
        for (int i=0; i<eigenValuesSqrt.length; i++) {
            for (int j=0; j<eigenValuesSqrt[i].length; j++) {
                eigenValuesSqrt[i][j] = Math.sqrt(newEigenValues[i][j]);
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
        normalDim = getNormalDim();

        latticePositions = new IVector[phase.getSpeciesMaster().moleculeCount()];

        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.space().makeVector();
            Atom atom = iterator.nextAtom();
            latticePositions[atomCount].E(atom.getType().getPositionDefinition().position(atom));
            atomCount++;
        }

        u = new double[normalDim];

        rRand = new double[waveVectors.length][normalDim];
        iRand = new double[waveVectors.length][normalDim];
        
        normalization = 1/Math.sqrt(phase.getSpeciesMaster().moleculeCount());
        
        // fills in elements of nominalU using NormalCoordWrapper
        iterator.reset();
        normalCoordMapper.setNumAtoms(iterator.size());
        iterator.reset();
        atomCount = 0;
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            normalCoordMapper.initNominalU(atom, atomCount);
            atomCount++;
        }
    }

    protected int getNormalDim() {
        //x, y, z
        // subclasses can override this to reserve space for other normal mode coordinates
        return phase.space().D();
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public boolean doTrial() {
        iterator.reset();
        int atomCount = 0;

        for (int iVector=0; iVector<waveVectors.length; iVector++) {
            for (int j=0; j<normalDim; j++) {
                rRand[iVector][j] = Simulation.random.nextGaussian() * eigenValuesSqrt[iVector][j];
                iRand[iVector][j] = Simulation.random.nextGaussian() * eigenValuesSqrt[iVector][j];
            }
        }
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            for (int i=0; i<normalDim; i++) {
                u[i] = 0;
            }
            for (int iVector=0; iVector<waveVectors.length; iVector++) {
                double kR = waveVectors[iVector].dot(latticePositions[atomCount]);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                

                for (int i=0; i<normalDim; i++) {
                    for (int j=0; j<normalDim; j++) {
                        u[j] += Math.sqrt(waveVectorCoefficients[iVector])*eigenVectors[iVector][i][j]*
                                  (rRand[iVector][i]*coskR - iRand[iVector][i]*sinkR);
                    }
                }
            }
            for (int i=0; i<normalDim; i++) {
                u[i] *= normalization;
            }
            normalCoordMapper.setToU(atom, atomCount, u);
            atomCount++;
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
    protected NormalCoordMapper normalCoordMapper;
    private final AtomIteratorAllMolecules iterator;
    private IVector[] latticePositions;
    private double[][] eigenValuesSqrt;
    private double[][][] eigenVectors;
    private IVector[] waveVectors;
    private double[] waveVectorCoefficients;
    protected int normalDim;
    protected double[] u;
    protected double[][] rRand;
    protected double[][] iRand;
    protected double normalization;
    public double msd;
}

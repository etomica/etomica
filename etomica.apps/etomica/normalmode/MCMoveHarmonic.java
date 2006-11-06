package etomica.normalmode;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.integrator.mcmove.MCMoveTracker;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;

public class MCMoveHarmonic extends MCMovePhase {

    public MCMoveHarmonic(PotentialMaster potentialMaster) {
        super(potentialMaster, new MCMoveTracker());
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setEigenValues(double[][] newEigenValues) {
        eigenValuesSqrt = new double[newEigenValues.length][newEigenValues[0].length];
        for (int i=0; i<eigenValuesSqrt.length; i++) {
            for (int j=0; j<eigenValuesSqrt[i].length; j++) {
                eigenValuesSqrt[i][j] = Math.sqrt(newEigenValues[i][j]);
            }
        }
    }
    
    public void setWaveVectors(Vector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }
    
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setPhase(Phase newPhase) {
        super.setPhase(newPhase);
        atomPos = phase.space().makeVector();
        atomActionTranslateTo = new AtomActionTranslateTo(phase.space());
        iterator.setPhase(newPhase);
        normalDim = getNormalDim();

        latticePositions = new Vector[phase.getSpeciesMaster().moleculeCount()];

        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.space().makeVector();
            Atom atom = iterator.nextAtom();
            latticePositions[atomCount].E(atom.type.getPositionDefinition().position(atom));
            atomCount++;
        }

        nominalU = new double[iterator.size()][normalDim];
        u = new double[normalDim];

        rRand = new double[waveVectors.length][normalDim];
        iRand = new double[waveVectors.length][normalDim];
        
        normalization = 2.0/Math.sqrt(phase.getSpeciesMaster().moleculeCount());
        
        // initialize what we think of as the original coordinates
        // allow sublcasses to initialize their own coordiantes
        initNominalU();
    }
    
    protected void initNominalU() {
        // fills in first D elements of nominalU with molecular x,y,z
        // subclasses can fill in other elements with their own
        // or not call this at all and not use x,y,z
        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            Vector moleculePos = atom.type.getPositionDefinition().position(atom);
            for (int i=0; i<moleculePos.D(); i++) {
                nominalU[atomCount][i] = moleculePos.x(i);
            }
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
                double g = Simulation.random.nextGaussian() * eigenValuesSqrt[iVector][j];
                double theta = Simulation.random.nextDouble()*2*Math.PI;
                rRand[iVector][j] = g*Math.cos(theta);
                iRand[iVector][j] = g*Math.sin(theta);
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
                        u[j] += eigenVectors[iVector][j][i]*(rRand[iVector][i]*coskR - iRand[iVector][i]*sinkR);
                    }
                }
            }
            setToU(atom, atomCount, normalization);
            atomCount++;
        }
        return true;
    }
    
    protected void setToU(Atom atom, int atomCount, double normalization) {
        for (int i=0; i<atomPos.D(); i++) {
            atomPos.setX(i, nominalU[atomCount][i] + u[i]*normalization);
        }
        atomActionTranslateTo.setDestination(atomPos);
        atomActionTranslateTo.actionPerformed(atom);
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
    private final AtomIteratorAllMolecules iterator;
    private AtomActionTranslateTo atomActionTranslateTo;
    private Vector[] latticePositions;
    private double[][] eigenValuesSqrt;
    private double[][][] eigenVectors;
    private Vector[] waveVectors;
    private Vector atomPos;
    protected int normalDim;
    protected double[][] nominalU;
    protected double[] u;
    protected double[][] rRand;
    protected double[][] iRand;
    protected double normalization;
    public double msd;
}

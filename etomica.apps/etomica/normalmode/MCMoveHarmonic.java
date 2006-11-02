package etomica.normalmode;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
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
    
    public void setEigenValues(double[] newEigenValues) {
        eigenValues = new double[newEigenValues.length];
        for (int i=0; i<eigenValues.length; i++) {
            eigenValues[i] = 0.5/newEigenValues[i];
        }
    }
    
    public void setWaveVectors(Vector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }
    
    public void setEigenVectors(double[][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setPhase(Phase newPhase) {
        super.setPhase(newPhase);
        iterator.setPhase(newPhase);
        normalDim = getNormalDim();

        latticePositions = new Vector[phase.getSpeciesMaster().moleculeCount()];

        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.space().makeVector();
            Atom atom = iterator.nextAtom();
            Vector atomPos = atom.type.getPositionDefinition().position(atom);
            latticePositions[atomCount].E(atomPos);
            atomCount++;
        }

        nominalU = new double[iterator.size()][normalDim];
        u = new double[normalDim];
        realT = new double[normalDim];
        imaginaryT = new double[normalDim];
        
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
            Vector atomPos = atom.type.getPositionDefinition().position(atom);
            for (int i=0; i<atomPos.D(); i++) {
                nominalU[atomCount][i] = atomPos.x(i);
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
        while (iterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
            Vector latticeSite = latticePositions[atom.getGlobalIndex()];
            atom.coord.position().E(latticeSite);
        }
        int D = potential.getSpace().D();
        Vector temp = potential.getSpace().makeVector();
        for (int i=0; i<eigenValues.length; i++) {
            double rand = Simulation.random.nextGaussian() * eigenValues[i];
            
            iterator.reset();
            int atomCount = 0;
            while (iterator.hasNext()) {
                temp.E(0);
                Atom atom = iterator.nextAtom();
                int offset = atomCount * D;
                for (int j=0; j<D; j++) {
                    temp.setX(j, eigenVectors[i][offset+j]);
                }
                atom.coord.position().PEa1Tv1(rand,temp);
            }
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
    private final AtomIteratorAllMolecules iterator;
    private Vector[] latticePositions;
    private double[] eigenValues;
    private double[][] eigenVectors;
    private Vector[] waveVectors;
    protected int normalDim;
    protected double[][] nominalU;
    protected double[] u;
}

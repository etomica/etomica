package etomica.normalmode;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.units.Energy;

/**
 * Meter that calculates the harmonic energy of a configuration given
 * eigenvectors and omegas corresponding to wave vectors.
 * @author Andrew Schultz
 */
public class MeterHarmonicEnergy extends DataSourceScalar {

    public MeterHarmonicEnergy() {
        super("Harmonic Energy", Energy.DIMENSION);
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    public double getDataAsScalar() {
        double energySum = 0;
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            for (int i=0; i<normalDim; i++) {
                realT[i] = 0;
                imaginaryT[i] = 0;
            }
            iterator.reset();
            int atomCount = 0;
            // sum T over atoms
            while (iterator.hasNext()) {
                Atom atom = iterator.nextAtom();
                coordinateDefinition.calcU(atom, atomCount, u);
                double kR = waveVectors[iVector].dot(latticePositions[atomCount]);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int i=0; i<normalDim; i++) {
                    realT[i] += coskR * u[i];
                    imaginaryT[i] += sinkR * u[i];
                }
                
                atomCount++;
            }
            
            // we want to calculate Q = A T
            // where A is made up of eigenvectors as columns
            for (int i=0; i<normalDim; i++) {
                double realCoord = 0, imaginaryCoord = 0;
                for (int j=0; j<normalDim; j++) {
                    realCoord += eigenVectors[iVector][j][i] * realT[j];
                    imaginaryCoord += eigenVectors[iVector][j][i] * imaginaryT[j];
                }
                // we were supposed to divide T by sqrt(atomCount), but it's easier to handle that here
                double normalCoord = (realCoord*realCoord + imaginaryCoord*imaginaryCoord)/atomCount;
                energySum += waveVectorCoefficients[iVector] * normalCoord * omegaSquared[iVector][i];
            }
        }
        return 0.5*energySum;
    }

    public Phase getPhase() {
        return phase;
    }

    public void setPhase(Phase newPhase) {
        phase = newPhase;
        iterator.setPhase(phase);
        normalDim = coordinateDefinition.getCoordinateDim();

        latticePositions = new IVector[phase.getSpeciesMaster().moleculeCount()];

        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.getSpace().makeVector();
            Atom atom = iterator.nextAtom();
            IVector atomPos = atom.getType().getPositionDefinition().position(atom);
            latticePositions[atomCount].E(atomPos);
            atomCount++;
        }

        coordinateDefinition.setNumAtoms(iterator.size());
        u = new double[normalDim];
        realT = new double[normalDim];
        imaginaryT = new double[normalDim];
        
        // notifies NormalCoordWrapper of the nominal position of each atom
        iterator.reset();
        atomCount = 0;
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            coordinateDefinition.initNominalU(atom, atomCount);
            atomCount++;
        }
    }
    
    public void setWaveVectors(IVector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    public void setEigenvectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] newOmegaSquared) {
        omegaSquared = newOmegaSquared;
    }
    
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    protected IVector[] latticePositions;
    protected final AtomIteratorAllMolecules iterator;
    protected Phase phase;
    protected int normalDim;
    protected double[] u;
    protected double[] realT, imaginaryT;
    protected IVector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenVectors;
    protected double[][] omegaSquared;
}

package etomica.normalmode;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.phase.Phase;
import etomica.space.Vector;
import etomica.units.Energy;

/**
 * Meter that calculates the harmonic energy of each normal mode for a 
 * configuration given eigenvectors and omegas corresponding to wave vectors.
 * @author Andrew Schultz
 */
public class MeterHarmonicSingleEnergy implements Meter {

    public MeterHarmonicSingleEnergy() {
        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{0});
        iterator = new AtomIteratorAllMolecules();
        tag = new DataTag();
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    

    public Data getData() {
        double[] x = data.getData();
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
                calcU(atom, atomCount);
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
                    realCoord += realT[i] * eigenVectors[iVector][j][i];
                    imaginaryCoord += imaginaryT[i] * eigenVectors[iVector][j][i];
                }
                double normalCoord = realCoord*realCoord + imaginaryCoord*imaginaryCoord;
                x[iVector*normalDim+i] = 0.5* normalCoord * omegaSquared[iVector][i];
            }
        }
        return data;
    }
    
    /**
     * Calculates the array of u elements for the given atom
     * subclasses should override this to fill in their own values
     */
    protected void calcU(Atom atom, int atomCount) {
        Vector pos = atom.type.getPositionDefinition().position(atom);
        for (int i=0; i<pos.D(); i++) {
            u[i] = pos.x(i) - nominalU[atomCount][i];
        }
    }

    public Phase getPhase() {
        return phase;
    }

    public void setPhase(Phase newPhase) {
        phase = newPhase;
        iterator.setPhase(phase);
        normalDim = getNormalDim();
        dataInfo = new DataInfoDoubleArray("Harmonic single energy", Energy.DIMENSION, new int[]{waveVectors.length,normalDim});
        data = new DataDoubleArray(new int[]{waveVectors.length,normalDim});

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
    
    public void setWaveVectors(Vector[] newWaveVectors) {
        waveVectors = newWaveVectors;
    }
    
    public void setEigenvectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] newOmegaSquared) {
        omegaSquared = newOmegaSquared;
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    private static final long serialVersionUID = 1L;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    private final DataTag tag;
    protected Vector[] latticePositions;
    protected double[][] nominalU;
    protected final AtomIteratorAllMolecules iterator;
    protected Phase phase;
    protected int normalDim;
    protected double[] u;
    protected double[] realT, imaginaryT;
    protected Vector[] waveVectors;
    protected double[][][] eigenVectors;
    protected double[][] omegaSquared;
    protected String name;
}

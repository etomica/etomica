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
    
    public void setNormalCoordMapper(NormalCoordMapper newNormalCoordWrapper) {
        normalCoordMapper = newNormalCoordWrapper;
        normalDim = normalCoordMapper.getNormalDim();
    }
    
    public NormalCoordMapper getNormalCoordWrapper() {
        return normalCoordMapper;
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
                normalCoordMapper.calcU(atom, atomCount, u);
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
                    realCoord += realT[j] * eigenVectors[iVector][j][i];
                    imaginaryCoord += imaginaryT[j] * eigenVectors[iVector][j][i];
                }
                // we were supposed to divide T by sqrt(atomCount), but it's easier to handle that here
                double normalCoord = (realCoord*realCoord + imaginaryCoord*imaginaryCoord)/atomCount;
                x[iVector*normalDim+i] = Math.exp(-0.5 * waveVectorCoefficients[iVector] * 
                        normalCoord * omegaSquared[iVector][i] / temperature);
            }
        }
        return data;
    }
    
    public Phase getPhase() {
        return phase;
    }

    public void setPhase(Phase newPhase) {
        phase = newPhase;
        iterator.setPhase(phase);
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

        normalCoordMapper.setNumAtoms(iterator.size());
        u = new double[normalDim];
        realT = new double[normalDim];
        imaginaryT = new double[normalDim];
        
        // fills in elements of nominalU using NormalCoordWrapper
        iterator.reset();
        atomCount = 0;
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            normalCoordMapper.initNominalU(atom, atomCount);
            atomCount++;
        }
    }
    
    public void setWaveVectors(Vector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    public void setEigenvectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] newOmegaSquared) {
        omegaSquared = newOmegaSquared;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    private static final long serialVersionUID = 1L;
    protected NormalCoordMapper normalCoordMapper;
    protected int normalDim;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    private final DataTag tag;
    protected Vector[] latticePositions;
    protected final AtomIteratorAllMolecules iterator;
    protected Phase phase;
    protected double temperature;
    protected double[] u;
    protected double[] realT, imaginaryT;
    protected Vector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected double[][][] eigenVectors;
    protected double[][] omegaSquared;
    protected String name;
}

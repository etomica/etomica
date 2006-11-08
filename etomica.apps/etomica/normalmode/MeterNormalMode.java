package etomica.normalmode;

import java.io.Serializable;

import etomica.action.Action;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.phase.Phase;
import etomica.space.Vector;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Null;

public class MeterNormalMode implements Meter, Action, Serializable {

    public MeterNormalMode() {
        tag = new DataTag();
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setNormalCoordWrapper(NormalCoordMapper newNormalCoordWrapper) {
        normalCoordMapper = newNormalCoordWrapper;
    }
    
    public NormalCoordMapper getNormalCoordWrapper() {
        return normalCoordMapper;
    }
    
    public void setWaveVectorFactory(WaveVectorFactory newWaveVectorFactory) {
        waveVectorFactory = newWaveVectorFactory;
    }
     
    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    /**
     * Sets the phase.  This method should be called when the Atoms are in 
     * their lattice positions.
     */
    public void setPhase(Phase newPhase) {
        callCount = 0;

        phase = newPhase;
        latticePositions = new Vector[phase.getSpeciesMaster().moleculeCount()];
        normalDim = normalCoordMapper.getNormalDim();

        iterator.setPhase(phase);
        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.space().makeVector();
            Atom atom = iterator.nextAtom();
            latticePositions[atomCount].E(atom.type.getPositionDefinition().position(atom));
            atomCount++;
        }

        waveVectorFactory.makeWaveVectors(phase);
        waveVectors = waveVectorFactory.getWaveVectors();
        // we don't actually care about the coefficients
        numWaveVectors = waveVectors.length;

        DataDoubleArray[] S = new DataDoubleArray[numWaveVectors];
        for (int i=0; i<S.length; i++) {
            // real and imaginary parts
            S[i] = new DataDoubleArray(new int[]{normalDim,normalDim});
        }
        data = new DataGroup(S);
        DataInfoDoubleArray[] Sinfo = new DataInfoDoubleArray[numWaveVectors];
        CompoundDimension area = new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2});
        for (int i=0; i<Sinfo.length; i++) {
            Sinfo[i] = new DataInfoDoubleArray("S", area, new int[]{normalDim,normalDim});
        }
        dataInfo = new DataInfoGroup("all S", Null.DIMENSION, Sinfo);

        u = new double[normalDim];
        realT = new double[normalDim];
        imaginaryT = new double[normalDim];
        
        // notifies NormalCoordWrapper of the nominal position of each atom
        normalCoordMapper.setNumAtoms(iterator.size());
        iterator.reset();
        atomCount = 0;
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            normalCoordMapper.initNominalU(atom, atomCount);
            atomCount++;
        }
    }
    
    public Phase getPhase() {
        return phase;
    }
    
    public Vector[] getWaveVectors() {
        return waveVectors;
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Calculating things and adds terms to the sums
     */
    public void actionPerformed() {
        callCount++;

        // |data.E(0)| here to calculate the current value rather than the sum
        // loop over wave vectors
        for (int iVector = 0; iVector < numWaveVectors; iVector++) {
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

            // add to S(k).  imaginary part of S is 0
            double[] sValues = ((DataDoubleArray)data.getData(iVector)).getData();
            for (int i=0; i<normalDim; i++) {
                for (int j=0; j<normalDim; j++) {
                    sValues[i*normalDim+j] += realT[i]*realT[j] + imaginaryT[i]*imaginaryT[j];
                }
            }
        }
    }

    /**
     * Returns the DataGroup of S(k) Tensors corresponding to the sum of 
     * T(k)*transpose(T(-k)).  To get the average (U), divide by 
     * (sqrt(numAtoms)*callCount())
     */
    public Data getData() {
        return data;
    }

    /**
     * Sets the tensor summation to 0.
     */
    public void reset() {
        data.E(0);
        callCount = 0;
    }
    
    public int getCallCount() {
        return callCount;
    }
    
    public void setName(String newName) {
        name = newName;
    }
    
    public String getName() {
        return name;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public String getLabel() {
        return "a label";
    }
    
    private static final long serialVersionUID = 1L;
    private Vector[] waveVectors;
    private WaveVectorFactory waveVectorFactory;
    protected NormalCoordMapper normalCoordMapper;
    private int numWaveVectors;
    private Phase phase;
    private String name;
    private final DataTag tag;
    private DataInfo dataInfo;
    private DataGroup data;
    private final AtomIteratorAllMolecules iterator;
    private Vector[] latticePositions;
    private int callCount;

    protected int normalDim;
    protected double[] u;
    protected double[] realT, imaginaryT;
}

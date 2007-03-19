package etomica.normalmode;

import java.io.Serializable;

import etomica.action.Action;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Null;

public class MeterNormalMode implements DataSource, Action, Serializable {

    public MeterNormalMode() {
        tag = new DataTag();
        iterator = new AtomIteratorAllMolecules();
    }
    
    /**
     * Sets the object that defines the real-space generalized coordinates.
     */
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
    }
    
    /**
     * @return the CoordinateDefinition last given via the set method.
     */
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }
    
    /**
     * Sets the object that defines the normal-coordinate wave vectors.
     */
    public void setWaveVectorFactory(WaveVectorFactory newWaveVectorFactory) {
        waveVectorFactory = newWaveVectorFactory;
    }

    /**
     * @return the WaveVectorFactory last given via the set methods.
     */
    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    /**
     * Sets the phase, and should be called while the Atoms are in 
     * their lattice positions.
     */
    public void setPhase(Phase newPhase) {
        callCount = 0;

        phase = newPhase;
        latticePositions = new IVector[phase.getSpeciesMaster().moleculeCount()];
        coordinateDim = coordinateDefinition.getCoordinateDim();

        iterator.setPhase(phase);
        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.getSpace().makeVector();
            Atom atom = iterator.nextAtom();
            latticePositions[atomCount].E(atom.getType().getPositionDefinition().position(atom));
            atomCount++;
        }

        waveVectorFactory.makeWaveVectors(phase);
        waveVectors = waveVectorFactory.getWaveVectors();
        // we don't actually care about the coefficients
        numWaveVectors = waveVectors.length;

        DataDoubleArray[] S = new DataDoubleArray[numWaveVectors];
        for (int i=0; i<S.length; i++) {
            // real and imaginary parts
            S[i] = new DataDoubleArray(new int[]{coordinateDim,coordinateDim});
        }
        data = new DataGroup(S);
        DataInfoDoubleArray[] Sinfo = new DataInfoDoubleArray[numWaveVectors];
        CompoundDimension area = new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2});
        for (int i=0; i<Sinfo.length; i++) {
            Sinfo[i] = new DataInfoDoubleArray("S", area, new int[]{coordinateDim,coordinateDim});
        }
        dataInfo = new DataInfoGroup("all S", Null.DIMENSION, Sinfo);

        u = new double[coordinateDim];
        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];
        
        // notifies CoordinateDefinition of the nominal position of each atom
        coordinateDefinition.setNumAtoms(iterator.size());
        iterator.reset();
        atomCount = 0;
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            coordinateDefinition.initNominalU(atom, atomCount);
            atomCount++;
        }
    }
    
    public Phase getPhase() {
        return phase;
    }
    
    public IVector[] getWaveVectors() {
        return waveVectors;
    }
    
    public IDataInfo getDataInfo() {
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
            for (int i=0; i<coordinateDim; i++) {
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
                for (int i=0; i<coordinateDim; i++) {
                    realT[i] += coskR * u[i];
                    imaginaryT[i] += sinkR * u[i];
                }
                atomCount++;
            }

            // add to S(k).  imaginary part of S is 0
            double[] sValues = ((DataDoubleArray)data.getData(iVector)).getData();
            for (int i=0; i<coordinateDim; i++) {
                for (int j=0; j<coordinateDim; j++) {
                    sValues[i*coordinateDim+j] += realT[i]*realT[j] + imaginaryT[i]*imaginaryT[j];
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
    private IVector[] waveVectors;
    private WaveVectorFactory waveVectorFactory;
    protected CoordinateDefinition coordinateDefinition;
    private int numWaveVectors;
    private Phase phase;
    private String name;
    private final DataTag tag;
    private IDataInfo dataInfo;
    private DataGroup data;
    private final AtomIteratorAllMolecules iterator;
    private IVector[] latticePositions;
    private int callCount;

    protected int coordinateDim;
    protected double[] u;
    protected double[] realT, imaginaryT;
}

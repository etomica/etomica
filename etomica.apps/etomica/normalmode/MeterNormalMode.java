package etomica.normalmode;

import etomica.action.Action;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.ConfigurationLattice;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Null;

public class MeterNormalMode implements Meter, Action {

    public MeterNormalMode() {
        tag = new DataTag();
        iterator = new AtomIteratorAllMolecules();
    }
    
    public void setPrimitive(Primitive newPrimitive) {
        primitive = newPrimitive;
    }
    
    public Primitive getPrimitive() {
        return primitive;
    }
    
    /**
     * Sets the phase.  This method should be called when the Atoms are in 
     * their lattice positions.
     */
    public void setPhase(Phase newPhase) {
        callCount = 0;

        phase = newPhase;
        latticePositions = new Vector[phase.getSpeciesMaster().moleculeCount()];
        normalDim = getNormalDim();

        iterator.setPhase(phase);
        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.space().makeVector();
            Atom atom = iterator.nextAtom();
            latticePositions[atomCount].E(atom.type.getPositionDefinition().position(atom));
            atomCount++;
        }

        if (primitive == null) {
            throw new RuntimeException("Please set primitive before the phase!!!!  Start again.");
        }
        
        numCells = 0;
        double d = -1;
        for (int i=0; i<phase.space().D(); i++) {
            //XXX divide by sqrt(2) for FCC
            int n = (int)Math.round(phase.getBoundary().getDimensions().x(i) / (primitive.getSize()[i]*Math.sqrt(2)));
            if (i>0 && n != numCells) {
                System.out.println("Things would be so much happier if you would just use the same number of cells in each direction.");
                throw new RuntimeException("... unless you're actually using a single-atom basis, in which case, you should fix the class.");
            }
            numCells = n;
            d = primitive.getSize()[i];
        }
        
        // FCC has 4-atom basis
        numWaveVectors = 4;
        for (int i=0; i<phase.space().D(); i++) {
            // -halfSize to halfSize in the other directions, including 0
            numWaveVectors *= numCells;
        }
//        numWaveVectors = (numWaveVectors-1)/2;
        
        //XXX the given constraints are for FCC
        Vector[] waveVectorsTemp = new Vector[numWaveVectors];
        int count = -1;
        for (int kx = 0; kx <= numCells; kx++) {
            for (int ky = ((kx==0) ? 1 : -numCells + 1); ky <= numCells; ky++) {
                for (int kz = ((kx==0 && ky==0) ? 1 : -numCells + 1); kz <= numCells; kz++) {
                    if (2 * (kx + ky + kz) <= 3 * numCells
                            && 2 * (kx + ky + kz) > -3 * numCells
                            && 2 * (kx + ky - kz) <= 3 * numCells
                            && 2 * (kx + ky - kz) > -3 * numCells
                            && 2 * (kx - ky + kz) <= 3 * numCells
                            && 2 * (kx - ky + kz) > -3 * numCells
                            && 2 * (kx - ky - kz) <= 3 * numCells
                            && 2 * (kx - ky - kz) > -3 * numCells) {
                        waveVectorsTemp[++count] = new Vector3D(kx, ky, kz);
                        waveVectorsTemp[count].TE(Math.sqrt(2) * Math.PI / d / numCells);
                    }
                }
            }
        }
        numWaveVectors = count+1;
        waveVectors = new Vector[numWaveVectors];
        System.arraycopy(waveVectorsTemp,0,waveVectors,0,numWaveVectors);

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

        nominalU = new double[iterator.size()][normalDim];
        u = new double[normalDim];
        realT = new double[normalDim];
        imaginaryT = new double[normalDim];
        
        // initialize what we think of as the original coordinates
        // allow sublcasses to initialize their own coordiantes
        initNominalU();
    }
    
    public Phase getPhase() {
        return phase;
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
//            if (iVector == 0) {
//                System.out.println("in mnm");
//                System.out.println("realT "+realT[0]+" "+realT[1]+" "+realT[2]);
//                System.out.println("imagT "+imaginaryT[0]+" "+imaginaryT[1]+" "+imaginaryT[2]);
//            }
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
     * Calculates the array of u elements for the given atom
     * subclasses should override this to fill in their own values
     */
    protected void calcU(Atom atom, int atomCount) {
        Vector pos = atom.type.getPositionDefinition().position(atom);
        for (int i=0; i<pos.D(); i++) {
            u[i] = pos.x(i) - nominalU[atomCount][i];
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
    
    private Vector[] waveVectors;
    private int numCells;
    private int numWaveVectors;
    private Phase phase;
    private String name;
    private final DataTag tag;
    private DataInfo dataInfo;
    private DataGroup data;
    private final AtomIteratorAllMolecules iterator;
    private Vector[] latticePositions;
    private int callCount;
    private Primitive primitive;

    protected double[][] nominalU;
    protected int normalDim;
    protected double[] u;
    protected double[] realT, imaginaryT;

    public static void main(String[] args) {
        MeterNormalMode foo = new MeterNormalMode();
        Simulation sim = new Simulation(Space3D.getInstance());
        Phase phase = new Phase(sim);
        Species species = new SpeciesSpheresMono(sim);
        phase.getAgent(species).setNMolecules(32);
        phase.setDimensions(Space.makeVector(new double[]{4, 4, 4}));
        LatticeCubicFcc lattice = new LatticeCubicFcc();
        ConfigurationLattice config = new ConfigurationLattice(lattice);
        foo.setPrimitive(lattice.getPrimitiveFcc());
        config.initializeCoordinates(phase);
        foo.setPhase(phase);
        foo.actionPerformed();
    }
}

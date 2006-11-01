package etomica.models.hexane;

import etomica.action.Action;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.ConfigurationLattice;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.meter.Meter;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTensor.DataInfoTensor;
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
        iterator = new AtomIteratorLeafAtoms();
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

        iterator.setPhase(phase);
        iterator.reset();
        int atomCount = 0;
        while (iterator.hasNext()) {
            latticePositions[atomCount] = phase.space().makeVector();
            latticePositions[atomCount].E(((AtomLeaf)iterator.nextAtom()).coord.position());
            atomCount++;
        }

        u = phase.space().makeVector();

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
        
        //XXX the given constraints are for FCC
        waveVectors = new Vector[numWaveVectors];
        int count = -1;
        for (int kx = -numCells + 1; kx <= numCells; kx++) {
            for (int ky = -numCells + 1; ky <= numCells; ky++) {
                for (int kz = -numCells + 1; kz <= numCells; kz++) {
                    if (2 * (kx + ky + kz) <= 3 * numCells
                            && 2 * (kx + ky + kz) > -3 * numCells
                            && 2 * (kx + ky - kz) <= 3 * numCells
                            && 2 * (kx + ky - kz) > -3 * numCells
                            && 2 * (kx - ky + kz) <= 3 * numCells
                            && 2 * (kx - ky + kz) > -3 * numCells
                            && 2 * (kx - ky - kz) <= 3 * numCells
                            && 2 * (kx - ky - kz) > -3 * numCells) {
                        waveVectors[++count] = new Vector3D(kx, ky, kz);
                        waveVectors[count].TE(Math.sqrt(2) * Math.PI / d / numCells);
                    }
                }
            }
        }

        realT = phase.space().makeVector();
        imaginaryT = phase.space().makeVector();
        DataTensor[] S = new DataTensor[numWaveVectors];
        for (int i=0; i<S.length; i++) {
            // real and imaginary parts
            S[i] = new DataTensor(phase.space());
        }
        data = new DataGroup(S);
        DataInfoTensor[] Sinfo = new DataInfoTensor[numWaveVectors];
        CompoundDimension area = new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2});
        for (int i=0; i<Sinfo.length; i++) {
            Sinfo[i] = new DataInfoTensor("S", area, phase.space());
        }
        dataInfo = new DataInfoGroup("all S", Null.DIMENSION, Sinfo);

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
            realT.E(0);
            imaginaryT.E(0);
            iterator.reset();
            int atomCount = 0;
            // sum T over atoms
            while (iterator.hasNext()) {
                AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
                Vector atomPos = atom.coord.position();
                u.E(atomPos);
                u.ME(latticePositions[atomCount]);
                double kR = waveVectors[iVector].dot(latticePositions[atomCount]);
                realT.PEa1Tv1(Math.cos(kR),u);
                imaginaryT.PEa1Tv1(Math.sin(kR),u);
                atomCount++;
            }
            // add to S(k).  imaginary part of S is 0
            DataTensor realS = (DataTensor)data.getData(iVector);
            realS.x.PEv1v2(realT, realT);
            realS.x.PEv1v2(imaginaryT, imaginaryT);
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
    private Vector realT, imaginaryT;
    private final AtomIteratorLeafAtoms iterator;
    private Vector[] latticePositions;
    private Vector u;
    private int callCount;
    private Primitive primitive;
    
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

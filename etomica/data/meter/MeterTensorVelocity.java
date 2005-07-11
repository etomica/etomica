package etomica.data.meter;
import etomica.Atom;
import etomica.AtomTypeLeaf;
import etomica.Data;
import etomica.DataInfo;
import etomica.EtomicaInfo;
import etomica.Meter;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceAtomic;
import etomica.data.types.DataTensor;
import etomica.space.ICoordinateKinetic;
import etomica.units.Dimension;

/**
 * A meter to compute the velocity component of the pressure tensor. 
 * Averages a tensor quantity formed from a dyad of the velocity of each atom. 
 * Specifically, the quantity averaged is 1/N * sum(pp/m), where p is the momentum,
 * m is the mass, and the sum is over all N atoms.
 * 
 * @author Rob Riggleman
 */

public class MeterTensorVelocity implements Meter, DataSourceAtomic, java.io.Serializable {
    /**
     * Iterator of atoms.
     */
    private final AtomIteratorPhaseDependent ai1 = new AtomIteratorLeafAtoms();
    
    public MeterTensorVelocity(Space space) {
        data = new DataTensor(space,"pp/m",Dimension.ENERGY);
        atomData = new DataTensor(space,"pp/m",Dimension.ENERGY);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Velocity tensor, formed from averaging dyad of velocity vector for each atom");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
       
    /**
     * Returns the dimension of the measured value, here given as energy
     */
    public Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Descriptive label
     *
     * @return "pp/m"
     */
    public String getLabel() {return "pp/m";}
    
    /**
     * Returns the velocity dyad (mass*vv) summed over all atoms, and divided by N
     */
    public Data getData() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        ai1.setPhase(phase);
        ai1.reset();
        data.E(0.0);
        int count = 0;
        while(ai1.hasNext()) {
            getData(ai1.nextAtom());
            data.PE(atomData);
            count++;
        }
        data.TE(1.0/count);
        return data;
    }
    
    /**
     * Returns the velocity dyad (mass*vv) for the given atom.
     */
    public Data getData(Atom atom) {
        atomData.x.E(((ICoordinateKinetic)atom.coord).velocity(), ((ICoordinateKinetic)atom.coord).velocity());
        atomData.TE(((AtomTypeLeaf)atom.type).rm());
        return atomData;
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    private String name;
    private Phase phase;
    private final DataTensor data, atomData;
}

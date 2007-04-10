package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.phase.Phase;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Temperature;

public class MeterPressureHardTensor implements DataSource, IntegratorHard.CollisionListener, java.io.Serializable {
    
    public MeterPressureHardTensor(Space space) {
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("PV/Nk",Temperature.DIMENSION, space);
        v = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure tensor measured via components of virial averaged over hard collisions");
        return info;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        if (phase == null || integratorHard == null) throw new IllegalStateException("must call setPhase and integrator before using meter");
        double t = integratorHard.getCurrentTime();
        data.x.TE(-1/((t-t0)*phase.getSpace().D()));
        t0 = t;

        //We're using the instantaneous velocity tensor with the average virial tensor
        //not quite right, but works out in the end.
        iterator.reset();
        while (iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            v.Ev1v2(((ICoordinateKinetic)a).getVelocity(), ((ICoordinateKinetic)a).getVelocity());
            v.TE((((AtomTypeLeaf)a.getType()).rm()));
            data.x.PE(v);
        }

        data.x.TE(1.0/phase.atomCount());
    
        return data;
    }
    
    public void collisionAction(IntegratorHard.Agent agent) {
        Tensor lcvt = agent.collisionPotential.lastCollisionVirialTensor();
        data.x.PE(lcvt);
    }
    
    public void setIntegrator(IntegratorHard newIntegrator) {
        if (newIntegrator == integratorHard) {
            return;
        }
        if (integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        if (newIntegrator == null) {
            phase = null;
            return;
        }
        phase = integratorHard.getPhase();
        integratorHard = newIntegrator;
        integratorHard.addCollisionListener(this);
    }
    
    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private static final long serialVersionUID = 1L;
    private double t0;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private Tensor v;
    private IntegratorHard integratorHard;
    private String name;
    private Phase phase;
    private final DataTensor data;
    private final IDataInfo dataInfo;
    protected final DataTag tag;
}

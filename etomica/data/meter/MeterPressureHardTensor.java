package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
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
        velocityTensor = space.makeTensor();
        v = space.makeTensor();
        timer = new DataSourceCountTime();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure tensor measured via components of virial averaged over hard collisions");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        if (phase == null || integratorHard == null) throw new IllegalStateException("must call setPhase and integrator before using meter");
        double t = timer.getDataAsScalar();
        if (t > t0) {
            //XXX Wrong!  you can't use the instantaneous velocity tensor with the average virial tensor
            velocityTensor.E(0.);
            iterator.reset();
            while (iterator.hasNext()) {
                AtomLeaf a = (AtomLeaf)iterator.nextAtom();
                v.Ev1v2(((ICoordinateKinetic)a.coord).velocity(), ((ICoordinateKinetic)a.coord).velocity());
                v.TE((((AtomTypeLeaf)a.getType()).rm()));
                velocityTensor.PE(v);
            }
            
            velocityTensor.TE(1.0/phase.atomCount());
                    
            data.x.TE(-1/((t-t0)*phase.space().D()*phase.atomCount()));
            data.x.PE(velocityTensor);
        
            t0 = t;
        }
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
            integratorHard.removeListener(timer);
        }
        timer.reset();
        if (newIntegrator == null) {
            phase = null;
            return;
        }
        phase = integratorHard.getPhase();
        integratorHard = newIntegrator;
        integratorHard.addCollisionListener(this);
        integratorHard.addListener(timer);
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
    private Tensor velocityTensor, v;
    private IntegratorHard integratorHard;
    private DataSourceCountTime timer;
    private String name;
    private Phase phase;
    private final DataTensor data;
    private final DataInfo dataInfo;
    protected final DataTag tag;
}

package etomica.data.meter;
import etomica.Atom;
import etomica.AtomTypeLeaf;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceCountTime;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorHard;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.units.Dimension;

public class MeterPressureHardTensor implements DataSource, IntegratorHard.CollisionListener {
    
    public MeterPressureHardTensor(Space space) {
        //XXX temperature, really?
        data = new DataTensor(space,new DataInfo("PV/Nk",Dimension.TEMPERATURE));
        velocityTensor = space.makeTensor();
        v = space.makeTensor();
        timer = new DataSourceCountTime();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure tensor measured via components of virial averaged over hard collisions");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public Data getData() {
        if (phase == null || integratorHard == null) throw new IllegalStateException("must call setPhase and integrator before using meter");
        double t = timer.getDataAsScalar();
        if (t > t0) {
            //XXX Wrong!  you can't use the instantaneous velocity tensor with the average virial tensor
            velocityTensor.E(0.);
            iterator.reset();
            while (iterator.hasNext()) {
                Atom a = iterator.nextAtom();
                v.E(((ICoordinateKinetic)a.coord).velocity(), ((ICoordinateKinetic)a.coord).velocity());
                v.TE((((AtomTypeLeaf)a.type).rm()));
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
    
    protected void setIntegrator(Integrator newIntegrator) {
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
        phase = integratorHard.getPhase()[0];
        integratorHard = (IntegratorHard)newIntegrator;
        integratorHard.addCollisionListener(this);
        integratorHard.addListener(timer);
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private double t0;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private Tensor velocityTensor, v;
    private IntegratorHard integratorHard;
    private DataSourceCountTime timer;
    private String name;
    private Phase phase;
    private DataTensor data;
}

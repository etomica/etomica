package etomica.modules.interfacial;
import etomica.api.IBox;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Temperature;

public class MeterVirialHardTensor implements DataSource, IntegratorHard.CollisionListener, java.io.Serializable {
    
    public MeterVirialHardTensor(ISpace space) {
    	dim = space.D();
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("PV/Nk",Temperature.DIMENSION, space);
        v = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        if (box == null || integratorHard == null) throw new IllegalStateException("must call setBox and integrator before using meter");
        double t = integratorHard.getCurrentTime();
        data.x.TE(-1/((t-t0)*dim));
        t0 = t;

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
        integratorHard = newIntegrator;
        if (newIntegrator == null) {
            box = null;
            return;
        }
        box = integratorHard.getBox();
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
    private Tensor v;
    private IntegratorHard integratorHard;
    private String name;
    private IBox box;
    private final DataTensor data;
    private final IDataInfo dataInfo;
    protected final DataTag tag;
    private final int dim;
}

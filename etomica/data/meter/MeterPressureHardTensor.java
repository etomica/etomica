package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.api.IAtomKinetic;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IData;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Temperature;

public class MeterPressureHardTensor implements IEtomicaDataSource, IntegratorHard.CollisionListener, java.io.Serializable {
    
    public MeterPressureHardTensor(ISpace space) {
    	dim = space.D();
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
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IData getData() {
        if (box == null || integratorHard == null) throw new IllegalStateException("must call setBox and integrator before using meter");
        double t = integratorHard.getCurrentTime();
        data.x.TE(-1/((t-t0)*dim));
        t0 = t;

        //We're using the instantaneous velocity tensor with the average virial tensor
        //not quite right, but works out in the end.
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            v.Ev1v2(a.getVelocity(), a.getVelocity());
            v.TE((((IAtomLeaf)a).getType().rm()));
            data.x.PE(v);
        }

        data.x.TE(1.0/box.getLeafList().getAtomCount());
    
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
    private final IEtomicaDataInfo dataInfo;
    protected final DataTag tag;
    private final int dim;
}

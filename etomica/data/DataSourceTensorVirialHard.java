package etomica.data;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Null;

/**
 * A MeterTensor that returns the virial component of the pressure tensor for a hard potential.  
 * This is the counterpart to MeterVelocityTensor, also called by the surface tension meter.
 *
 * @author Rob Riggleman
 */

public class DataSourceTensorVirialHard implements DataSource, EtomicaElement, IntegratorHard.CollisionListener, java.io.Serializable {
    
    public DataSourceTensorVirialHard(Space space) {
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("PV/NkT", Null.DIMENSION, space);
        timer = new DataSourceCountTime();
        work = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public DataSourceTensorVirialHard(Space space, IntegratorHard integrator) {
        this(space);
        setIntegrator(integrator);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Virial tensor for a hard potential");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    //TODO consider how to ensure timer is advanced before this method is invoked
    public Data getData() {
        double elapsedTime = timer.getDataAsScalar();
        if(elapsedTime == 0.0) {
            data.E(Double.NaN);
            return data;
        }
        Phase phase = integratorHard.getPhase();
        int D = phase.space().D();

        work.TE(-1./(integratorHard.getTemperature()*elapsedTime*D*phase.atomCount()));
        data.x.E(work);
        //don't add 1.0 to diagonal elements because meter returns only virial contribution to pressure
        timer.reset();
        return data;
    }
    
    /**
     * Sums contribution to virial for each collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        work.PE(agent.collisionPotential.lastCollisionVirialTensor());
    }
    
    /**
     * Contribution to the virial from the most recent collision of the given pair/potential.
     */
    public Tensor collisionValue(IntegratorHard.Agent agent) {
        data.x.E(agent.collisionPotential.lastCollisionVirialTensor());
        data.x.TE(1/(double)integratorHard.getPhase().atomCount());
        return data.x;
    }
                
    /**
     * Informs meter of the integrator for its phase, and zeros elapsed-time counter
     */
    public void setIntegrator(IntegratorHard newIntegrator) {
        if(newIntegrator == integratorHard) return;
        if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
            integratorHard.removeListener(timer);
        }
        integratorHard = newIntegrator;
        if(newIntegrator != null) {
            timer.reset();
            newIntegrator.addCollisionListener(this);
            newIntegrator.addListener(timer);
        }
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
    
    private String name;
    private DataSourceCountTime timer;
    private IntegratorHard integratorHard;
    private final DataTensor data;
    private final Tensor work;
    private final DataInfoTensor dataInfo;
    protected final DataTag tag;
}

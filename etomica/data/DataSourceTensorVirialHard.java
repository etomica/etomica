package etomica.data;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.integrator.IntegratorHard;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Null;

/**
 * A MeterTensor that returns the virial component of the pressure tensor for a hard potential.  
 * This is the counterpart to MeterVelocityTensor, also called by the surface tension meter.
 *
 * @author Rob Riggleman
 */

public class DataSourceTensorVirialHard implements DataSource, IntegratorHard.CollisionListener, java.io.Serializable {
    
    public DataSourceTensorVirialHard(Space space) {
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("PV/NkT", Null.DIMENSION, space);
        work = space.makeTensor();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public DataSourceTensorVirialHard(Space space, IntegratorHard integrator) {
        this(space);
        setIntegrator(integrator);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    public Data getData() {
        double currentTime = integratorHard.getCurrentTime();
        double elapsedTime = currentTime - lastTime;
        if(elapsedTime == 0.0) {
            data.E(Double.NaN);
            return data;
        }
        Box box = integratorHard.getBox();
        int D = box.getSpace().D();

        work.TE(-1./(integratorHard.getTemperature()*elapsedTime*D*box.atomCount()));
        data.x.E(work);
        //don't add 1.0 to diagonal elements because meter returns only virial contribution to pressure
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
        data.x.TE(1/(double)integratorHard.getBox().atomCount());
        return data.x;
    }
                
    /**
     * Informs meter of the integrator for its box, and zeros elapsed-time counter
     */
    public void setIntegrator(IntegratorHard newIntegrator) {
        if(newIntegrator == integratorHard) return;
        if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        integratorHard = newIntegrator;
        if(newIntegrator != null) {
            newIntegrator.addCollisionListener(this);
            lastTime = newIntegrator.getCurrentTime();
        }
    }
    
    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
    
    private static final long serialVersionUID = 1L;
    protected double lastTime;
    protected IntegratorHard integratorHard;
    protected final DataTensor data;
    protected final Tensor work;
    protected final DataInfoTensor dataInfo;
    protected final DataTag tag;
}

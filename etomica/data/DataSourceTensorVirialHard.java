package etomica.data;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Space;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorHard;
import etomica.space.Tensor;
import etomica.units.Dimension;

/**
 * A MeterTensor that returns the virial component of the pressure tensor for a hard potential.  
 * This is the counterpart to MeterVelocityTensor, also called by the surface tension meter.
 *
 * @author Rob Riggleman
 */

public class DataSourceTensorVirialHard implements DataSource, EtomicaElement, IntegratorHard.CollisionListener, java.io.Serializable {
    
    public DataSourceTensorVirialHard(Space space, IntegratorHard integrator) {
        data = new DataTensor(space,new DataInfo("PV/NkT",Dimension.NULL));
        lastData = new DataTensor(space,new DataInfo("PV/NkT",Dimension.NULL));
        timer = new DataSourceCountTime();
        setIntegrator(integrator);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Virial tensor for a hard potential");
        return info;
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    /**
     * Current value of the meter, obtained by dividing sum of collision virial contributions by time elapsed since last call.
     * If elapsed-time interval is zero, returns the value reported at the last call to the method.
     */
    //TODO consider how to ensure timer is advanced before this method is invoked
    public Data getData() {
        double elapsedTime = timer.getDataAsScalar();
        if(elapsedTime == 0.0) {
            lastData.E(Double.NaN);
            return lastData;
        }
        int D = phase.space().D();

        lastData.E(data);
        lastData.TE(-1./(integratorHard.getTemperature()*elapsedTime*D*phase.atomCount()));
        data.E(0);
        //don't add 1.0 to diagonal elements because meter returns only virial contribution to pressure
        timer.reset();
        return lastData;
    }
    
    /**
     * Sums contribution to virial for each collision.
     */
    public void collisionAction(IntegratorHard.Agent agent) {
       lastData.x.E(agent.collisionPotential.lastCollisionVirialTensor());
       data.PE(lastData);
    }
    
    /**
     * Contribution to the virial from the most recent collision of the given pair/potential.
     */
    public Tensor collisionValue(IntegratorHard.Agent agent) {
        lastData.x.E(agent.collisionPotential.lastCollisionVirialTensor());
        lastData.x.TE(1/(double)phase.atomCount());
        return lastData.x;
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
            integratorHard.addCollisionListener(this);
            newIntegrator.addListener(timer);
            setPhase(integratorHard.getPhase()[0]);
        }
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
    private void setPhase(Phase phase) {
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
    private DataSourceCountTime timer;
    private IntegratorHard integratorHard;
    private final DataTensor data, lastData;
}//end of MeterTensorVirialHard

package etomica.data.meter;

import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.phase.Phase;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterKineticEnergyFromIntegrator extends DataSourceScalar implements Meter, IntegratorNonintervalListener, java.io.Serializable {

    public MeterKineticEnergyFromIntegrator(IntegratorMD aIntegrator) {
        super("Kinetic Energy",Dimension.ENERGY);
        integrator = aIntegrator;
        integrator.addListener(this);
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public double getDataAsScalar() {
        return integrator.getKineticEnergy();
    }
 
    public int getPriority() {return 200;}
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorEvent.DONE) {
            double currentKE = integrator.getKineticEnergy();
            MeterKineticEnergy meterKE = new MeterKineticEnergy();
            meterKE.setPhase(phase);
            double KE = meterKE.getDataAsScalar();
            if (Math.abs(KE - currentKE) > 1.e-9*Math.abs(KE+currentKE)) {
                System.out.println("final kinetic energy ("+currentKE+") for "+phase+" doesn't match actual energy ("+KE+")");
                meterKE.getData();
            }
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
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private IntegratorMD integrator;
}

package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.phase.Phase;
import etomica.units.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterKineticEnergyFromIntegrator extends DataSourceScalar implements IntegratorNonintervalListener, java.io.Serializable {

    public MeterKineticEnergyFromIntegrator() {
        super("Kinetic Energy",Energy.DIMENSION);
    }
    
    public MeterKineticEnergyFromIntegrator(IntegratorMD aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }
    
    public void setIntegrator(IntegratorMD newIntegrator) {
        if (integrator != null) {
            integrator.removeListener(this);
        }
        integrator = newIntegrator;
        if (integrator != null) {
            integrator.addListener(this);
        }
    }
    
    public IntegratorMD getIntegrator() {
        return integrator;
    }
    
    public double getDataAsScalar() {
        return integrator.getKineticEnergy();
    }
 
    public int getPriority() {return 200;}
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorNonintervalEvent.DONE) {
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

    private static final long serialVersionUID = 1L;
    private Phase phase;
    private IntegratorMD integrator;
}

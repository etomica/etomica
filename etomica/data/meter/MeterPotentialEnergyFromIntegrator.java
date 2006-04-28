package etomica.data.meter;

import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.units.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterPotentialEnergyFromIntegrator extends DataSourceScalar implements Meter, IntegratorNonintervalListener, java.io.Serializable {

    public MeterPotentialEnergyFromIntegrator() {
        super("Potential Energy",Energy.DIMENSION);
    }
    
    public MeterPotentialEnergyFromIntegrator(IntegratorPhase aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }
    
    public void setIntegrator(IntegratorPhase newIntegrator) {
        if (integrator != null) {
            integrator.removeListener(this);
        }
        integrator = newIntegrator;
        if (integrator != null) {
            integrator.addListener(this);
        }
    }
    
    public IntegratorPhase getIntegrator() {
        return integrator;
    }
    
    public double getDataAsScalar() {
        return integrator.getPotentialEnergy();
    }
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorNonintervalEvent.DONE) {
            double currentPE = integrator.getPotentialEnergy();
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(integrator.getPotential());
            meterPE.setPhase(integrator.getPhase());
            double PE = meterPE.getDataAsScalar();
            if (Math.abs(PE - currentPE) > 1.e-9*Math.abs(PE+currentPE)) {
                System.out.println("final potential energy ("+currentPE+") for "+integrator.getPhase()+" doesn't match actual energy ("+PE+")");
                meterPE.getData();
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
    private IntegratorPhase integrator;
}

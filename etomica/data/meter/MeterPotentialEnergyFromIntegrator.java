package etomica.data.meter;

import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterPotentialEnergyFromIntegrator extends DataSourceScalar implements Meter, IntegratorNonintervalListener, java.io.Serializable {

    public MeterPotentialEnergyFromIntegrator(IntegratorPhase aIntegrator) {
        super("Potential Energy",Dimension.ENERGY);
        integrator = aIntegrator;
        integrator.addListener(this);
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public double getDataAsScalar() {
        return integrator.getPotentialEnergy();
    }
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorEvent.DONE) {
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

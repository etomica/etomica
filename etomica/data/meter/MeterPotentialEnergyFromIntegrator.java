package etomica.data.meter;

import etomica.Phase;
import etomica.data.DataInfo;
import etomica.data.DataSourceScalar;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterPotentialEnergyFromIntegrator extends DataSourceScalar implements Meter, IntegratorNonintervalListener, java.io.Serializable {

    public MeterPotentialEnergyFromIntegrator(Integrator aIntegrator) {
        super("Potential Energy",Dimension.ENERGY);
        integrator = aIntegrator;
        integrator.addListener(this);
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public double getDataAsScalar() {
        final double[] PE = integrator.getPotentialEnergy();
        final Phase[] integratorPhases = integrator.getPhase();
        for (int i=0; i<PE.length; i++) {
            if (integratorPhases[i] == phase) {
                return PE[i];
            }
        }
        throw new IllegalStateException("Meter's phase not handled by the Meter's integrator");
    }
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorEvent.DONE) {
            double[] currentPE = integrator.getPotentialEnergy();
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(integrator.getPotential());
            meterPE.setPhase(phase);
            double PE = meterPE.getDataAsScalar();
            final Phase[] integratorPhases = integrator.getPhase();
            for (int i=0; i<currentPE.length; i++) {
                if (integratorPhases[i] == phase) {
                    if (Math.abs(PE - currentPE[i]) > 1.e-9*Math.abs(PE+currentPE[i])) {
                        System.out.println("final potential energy ("+currentPE[i]+") for "+integratorPhases[i]+" doesn't match actual energy ("+PE+")");
                        meterPE.getData();
                    }
                    return;
                }
            }
            throw new IllegalStateException("Meter's phase not handled by the Meter's integrator");
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
    private Integrator integrator;
}

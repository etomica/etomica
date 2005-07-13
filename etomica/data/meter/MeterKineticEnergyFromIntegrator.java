package etomica.data.meter;

import etomica.DataInfo;
import etomica.IntegratorEvent;
import etomica.IntegratorNonintervalEvent;
import etomica.IntegratorNonintervalListener;
import etomica.Meter;
import etomica.Phase;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
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
        final double[] KE = integrator.getKineticEnergy();
        final Phase[] integratorPhases = integrator.getPhase();
        for (int i=0; i<KE.length; i++) {
            if (integratorPhases[i] == phase) {
                return KE[i];
            }
        }
        throw new IllegalStateException("Meter's phase not handled by the Meter's integrator");
    }
 
    public int getPriority() {return 200;}
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorEvent.DONE) {
            double[] currentKE = integrator.getKineticEnergy();
            MeterKineticEnergy meterKE = new MeterKineticEnergy();
            meterKE.setPhase(phase);
            double KE = meterKE.getDataAsScalar();
            final Phase[] integratorPhases = integrator.getPhase();
            for (int i=0; i<currentKE.length; i++) {
                if (integratorPhases[i] == phase) {
                    if (Math.abs(KE - currentKE[i]) > 1.e-9*Math.abs(KE+currentKE[i])) {
                        System.out.println("final kinetic energy ("+currentKE[i]+") for "+integratorPhases[i]+" doesn't match actual energy ("+KE+")");
                        meterKE.getData();
                    }
                    return;
                }
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

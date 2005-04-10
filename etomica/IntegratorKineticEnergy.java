package etomica;

import etomica.Integrator.IntervalEvent;
import etomica.data.meter.MeterKineticEnergy;
import etomica.integrator.IntegratorMD;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class IntegratorKineticEnergy implements DataSource, Integrator.IntervalListener {

    public IntegratorKineticEnergy(IntegratorMD aIntegrator) {
        integrator = aIntegrator;
        integrator.addIntervalListener(this);
        meterKE = new MeterKineticEnergy();
    }
    
    public double[] getData() {
        return integrator.getKineticEnergy();
    }
 
    /**
     * Length of data is the number of phases tracked by integrator; returns it.
     */
    public int getDataLength() {
        return integrator.getPhase().length;
    }

    public int getPriority() {return 200;}
    
    public void intervalAction(IntervalEvent evt) {
        if (evt.type() == Integrator.IntervalEvent.DONE) {
            double[] currentKE = integrator.getKineticEnergy();
            Phase[] phase = integrator.getPhase();
            meterKE.setPhase(phase);
            double[] KE = meterKE.getData();
            for (int i=0; i<phase.length; i++) {
                if (Math.abs(KE[i] - currentKE[i]) > 1.e-9) {
                    System.out.println("final kinetic energy ("+currentKE[i]+") for "+phase[i]+" doesn't match actual energy ("+KE[i]+")");
                }
            }
        }
    }
    
    public String getLabel() {
        return label;
    }
    public void setLabel(String string) {
        label = string;
    }

    public Dimension getDimension() {
        return Dimension.ENERGY;
    }

    public DataTranslator getTranslator() {
        return DataTranslator.IDENTITY;
    }

    private IntegratorMD integrator;
    private String label = "Kinetic Energy";
    private MeterKineticEnergy meterKE;
}

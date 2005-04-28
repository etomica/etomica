package etomica;

import etomica.data.meter.MeterPotentialEnergy;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class IntegratorPotentialEnergy implements DataSource, IntegratorNonintervalListener {

    public IntegratorPotentialEnergy(Integrator aIntegrator) {
        integrator = aIntegrator;
        integrator.addListener(this);
    }
    
    public double[] getData() {
        return integrator.getPotentialEnergy();
    }
    
    /**
     * Length of data is the number of phases tracked by integrator; returns it.
     */
    public int getDataLength() {
        return integrator.getPhase().length;
    }

    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorIntervalEvent.DONE) {
            double[] currentPE = integrator.getPotentialEnergy();
            Phase[] phase = integrator.getPhase();
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(integrator.potential);
            meterPE.setPhase(phase);
            double[] PE = meterPE.getData();
            for (int i=0; i<phase.length; i++) {
                if (Math.abs(PE[i] - currentPE[i]) > 1.e-9*Math.abs(PE[i]+currentPE[i])) {
                    System.out.println("final potential energy ("+currentPE[i]+") for "+phase[i]+" doesn't match actual energy ("+PE[i]+")");
                    meterPE.getData();
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

    private Integrator integrator;
    private String label = "Potential Energy";
}

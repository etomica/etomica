package etomica;

import etomica.Integrator.IntervalEvent;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class IntegratorPotentialEnergy implements DataSource, Integrator.IntervalListener {

    public IntegratorPotentialEnergy(Integrator aIntegrator) {
        integrator = aIntegrator;
        integrator.addIntervalListener(this);
    }
    
    public double[] getData() {
        return integrator.getPotentialEnergy();
    }

    public int getPriority() {return 200;}
    
    public void intervalAction(IntervalEvent evt) {
        if (evt.type() == Integrator.IntervalEvent.DONE) {
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

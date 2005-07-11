package etomica;

import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class IntegratorPotentialEnergy implements DataSource, IntegratorNonintervalListener, java.io.Serializable {

    public IntegratorPotentialEnergy(Integrator aIntegrator) {
        data = new DataDoubleArray("Potential Energy",Dimension.ENERGY);
        integrator = aIntegrator;
        integrator.addListener(this);
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public Data getData() {
        final double[] x = data.getData();
        final double[] PE = integrator.getPotentialEnergy();
        for (int i=0; i<x.length; i++) {
            x[i] = PE[i];
        }
        return data;
    }
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorIntervalEvent.START) {
            data.setLength(evt.getSource().getPhase().length);
        }
        else if (evt.type() == IntegratorIntervalEvent.DONE) {
            double[] currentPE = integrator.getPotentialEnergy();
            Phase[] phase = integrator.getPhase();
            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(integrator.potential);
            for (int i=0; i<phase.length; i++) {
                meterPE.setPhase(phase[i]);
                double PE = meterPE.getDataAsScalar();
                if (Math.abs(PE - currentPE[i]) > 1.e-9*Math.abs(PE+currentPE[i])) {
                    System.out.println("final potential energy ("+currentPE[i]+") for "+phase[i]+" doesn't match actual energy ("+PE+")");
                    meterPE.getData();
                }
            }
        }
    }

    private final DataDoubleArray data;
    private Integrator integrator;
}

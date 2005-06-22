package etomica;

import etomica.data.meter.MeterKineticEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.IntegratorMD;
import etomica.units.Dimension;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class IntegratorKineticEnergy implements DataSource, IntegratorNonintervalListener {

    public IntegratorKineticEnergy(IntegratorMD aIntegrator) {
        data = new DataDoubleArray(new DataInfo("Potential Energy",Dimension.ENERGY));
        integrator = aIntegrator;
        integrator.addListener(this);
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public Data getData() {
        final double[] x = data.getData();
        final double[] PE = integrator.getKineticEnergy();
        for (int i=0; i<x.length; i++) {
            x[i] = PE[i];
        }
        return data;
    }
 
    public int getPriority() {return 200;}
    
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorEvent.START) {
            data.setLength(evt.getSource().getPhase().length);
        }
        else if (evt.type() == IntegratorEvent.DONE) {
            double[] currentKE = integrator.getKineticEnergy();
            Phase[] phase = integrator.getPhase();
            MeterKineticEnergy meterKE = new MeterKineticEnergy();
            for (int i=0; i<phase.length; i++) {
                meterKE.setPhase(phase[i]);
                double KE = meterKE.getDataAsScalar();
                if (Math.abs(KE - currentKE[i]) > 1.e-9) {
                    System.out.println("final kinetic energy ("+currentKE[i]+") for "+phase[i]+" doesn't match actual energy ("+KE+")");
                }
            }
        }
    }
    
    private final DataDoubleArray data;
    private IntegratorMD integrator;
}

package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
import etomica.units.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterPotentialEnergyFromIntegrator extends DataSourceScalar {

    public MeterPotentialEnergyFromIntegrator() {
        super("Potential Energy",Energy.DIMENSION);
    }
    
    public MeterPotentialEnergyFromIntegrator(IntegratorBox aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }
    
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
    }
    
    public IntegratorBox getIntegrator() {
        return integrator;
    }
    
    public double getDataAsScalar() {
        return integrator.getPotentialEnergy();
    }
    
    private static final long serialVersionUID = 1L;
    private IntegratorBox integrator;
}

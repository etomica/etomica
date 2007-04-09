package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorPhase;
import etomica.units.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterPotentialEnergyFromIntegrator extends DataSourceScalar {

    public MeterPotentialEnergyFromIntegrator() {
        super("Potential Energy",Energy.DIMENSION);
    }
    
    public MeterPotentialEnergyFromIntegrator(IntegratorPhase aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }
    
    public void setIntegrator(IntegratorPhase newIntegrator) {
        integrator = newIntegrator;
    }
    
    public IntegratorPhase getIntegrator() {
        return integrator;
    }
    
    public double getDataAsScalar() {
        return integrator.getPotentialEnergy();
    }
    
    private static final long serialVersionUID = 1L;
    private IntegratorPhase integrator;
}

package etomica.normalmode;

import etomica.api.IBox;
import etomica.data.DataSourceScalar;
import etomica.units.Energy;

/**
 * Meter that returns the potential energy using an IAPIPotential.
 * 
 * @author Andrew Schultz
 */
public class MeterAPIPotentialEnergy extends DataSourceScalar {
    
    public MeterAPIPotentialEnergy(IAPIPotential potential) {
        super("energy", Energy.DIMENSION);
        this.potential = potential;
    }

    public double getDataAsScalar() {
        return potential.calculateEnergy(box);
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }
    
    public IBox getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    protected final IAPIPotential potential;
    protected IBox box;
}

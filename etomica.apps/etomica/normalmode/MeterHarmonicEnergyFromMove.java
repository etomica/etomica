package etomica.normalmode;

import etomica.data.DataSourceScalar;
import etomica.units.Energy;

/**
 * Meter that returns the harmonic energy from the MC move as a data source.
 * 
 * @author Andrew Schultz
 */
public class MeterHarmonicEnergyFromMove extends DataSourceScalar {
    
    private static final long serialVersionUID = 1L;

    public MeterHarmonicEnergyFromMove(MCMoveHarmonic mcMoveHarmonic) {
        super("energy", Energy.DIMENSION);
        this.mcMoveHarmonic = mcMoveHarmonic;
    }

    public double getDataAsScalar() {
        return mcMoveHarmonic.getLastTotalEnergy();
    }

    protected final MCMoveHarmonic mcMoveHarmonic;
}

package etomica.modules.materialfracture;
import etomica.api.IBox;
import etomica.data.DataSourceScalar;
import etomica.units.Pressure2D;

/**
 * Meter that calculates the stress within the gage cell.
 *
 * @author Andrew Schultz
 */
public final class MeterLoad extends DataSourceScalar {
    
    public MeterLoad(PotentialCalculationForceStress pc) {
        super("Load", Pressure2D.DIMENSION);
        this.pc = pc;
    }

    public double getDataAsScalar(){
        return pc.getLoad();
    }

    protected final PotentialCalculationForceStress pc;
    protected IBox box;
}
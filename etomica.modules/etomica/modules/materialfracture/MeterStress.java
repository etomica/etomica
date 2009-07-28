package etomica.modules.materialfracture;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.units.Pressure2D;

/**
 * Meter that calculates the stress within the gage cell.
 *
 * @author Andrew Schultz
 */
public class MeterStress extends DataSourceScalar {
    
    public MeterStress(PotentialCalculationForceStress pc) {
        super("Stress", Pressure2D.DIMENSION);
        this.pc = pc;
    }

    public void setBox(IBox newBox) {
        box = newBox;
    }

    public IBox getBox() {
        return box;
    }

    public double getDataAsScalar(){
        double area = 1;
        IVector dim = box.getBoundary().getDimensions();
        for (int i=1; i<dim.getD(); i++) {
            area *= dim.getX(i);
        }

        return pc.getLoad()/area / 2.0;
    }

    protected final PotentialCalculationForceStress pc;
    protected IBox box;
}
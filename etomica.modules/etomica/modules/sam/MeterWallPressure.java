package etomica.modules.sam;

import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.units.Pressure;

public class MeterWallPressure extends DataSourceScalar {


    public MeterWallPressure(PotentialCalculationForceSumWall pc) {
        super("Wall Pressure", Pressure.DIMENSION);
        this.pc = pc;
    }

    public void setBox(IBox newBox) {
        box = newBox;
    }
    
    public double getDataAsScalar() {
        double f = pc.getWallForce();
        IVector dimensions = box.getBoundary().getDimensions();
        double A = 1;
        for (int i=0; i<dimensions.getD(); i++) {
            if (i == pc.getWallPotential().getWallDim()) {
                continue;
            }
            A *= dimensions.x(i);
        }
        return f / A;
    }

    private static final long serialVersionUID = 1L;
    protected final PotentialCalculationForceSumWall pc;
    protected IBox box;
}

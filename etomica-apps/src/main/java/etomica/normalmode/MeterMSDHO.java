package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.units.dimensions.Null;

public class MeterMSDHO extends DataSourceScalar {
    protected Box box;
    protected int nBeads;

    public MeterMSDHO(int nBeads, Box box) {
        super("Stuff", Null.DIMENSION);
        this.box = box;
        this.nBeads = nBeads;
    }

    @Override
    public double getDataAsScalar() {
        double r2sum = 0;
        for (int i=0; i < nBeads; i++){
            r2sum += box.getLeafList().get(i).getPosition().squared();
        }
        return r2sum/nBeads;
    }
}

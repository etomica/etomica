package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.units.dimensions.Null;

public class MeterPrimPI extends DataSourceScalar {
    protected Box box;
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN;

    public MeterPrimPI(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN) {
        super("Stuff", Null.DIMENSION);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
    }

    @Override
    public double getDataAsScalar() {

        pmBonding.computeAll(false);
        pcP1.computeAll(false);

        return 1.0/2.0/betaN + pcP1.getLastEnergy()- pmBonding.getLastEnergy();
    }
}

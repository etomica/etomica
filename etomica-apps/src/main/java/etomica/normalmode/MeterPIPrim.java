package etomica.normalmode;

import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.units.dimensions.Null;

public class MeterPIPrim extends DataSourceScalar {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;

    public MeterPIPrim(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads) {
        super("Stuff", Null.DIMENSION);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
    }

    @Override
    public double getDataAsScalar() {
        pmBonding.computeAll(false);
        pcP1.computeAll(false);
        double En_prim = 1.0/2.0/betaN + pcP1.getLastEnergy() - pmBonding.getLastEnergy();
//        System.out.println("prim: " + En_prim);
        return En_prim;
    }
}

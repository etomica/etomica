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
    protected int nBeads;

    public MeterPrimPI(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads) {
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
//        System.out.println(1.0/2.0/betaN + " " + pcP1.getLastEnergy()/nBeads + " " + pmBonding.getLastEnergy());
//        System.out.println((1.0/2.0/betaN - pmBonding.getLastEnergy()) + "  " + pcP1.getLastEnergy()/nBeads
//        + "    " + (1.0/2.0/betaN + pcP1.getLastEnergy()/nBeads - pmBonding.getLastEnergy()));
        return 1.0/2.0/betaN + pcP1.getLastEnergy()/nBeads - pmBonding.getLastEnergy();
    }
}

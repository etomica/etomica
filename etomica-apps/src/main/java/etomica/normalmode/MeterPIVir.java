package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIVir extends DataSourceScalar {
    protected Box box;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;

    public MeterPIVir(PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
        super("Stuff", Null.DIMENSION);
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;
    }

    @Override
    public double getDataAsScalar() {
        pcP1.computeAll(true);
        Vector[] forces = pcP1.getForces();
        double vir = 0;
        for (int i = 0; i < nBeads; i++){
            Vector xi = box.getLeafList().get(i).getPosition();
            vir -= forces[i].dot(xi);
        }
        double En_vir = pcP1.getLastEnergy()/nBeads + 1.0/2.0/nBeads*vir;
//        System.out.println("vir: " + En_vir);
        return En_vir;
    }
}

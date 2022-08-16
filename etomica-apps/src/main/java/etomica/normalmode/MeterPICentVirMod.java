package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPICentVirMod extends DataSourceScalar {
    protected Box box;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;

    public MeterPICentVirMod(PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
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
        Vector xc = box.getSpace().makeVector();
        for (int i = 0; i < nBeads; i++){
            Vector xi = box.getLeafList().get(i).getPosition();
            xc.PE(xi);
        }
        xc.TE(1.0/nBeads);

        for (int i = 0; i < nBeads; i++){
            Vector xi = box.getLeafList().get(i).getPosition();
            vir -= forces[i].dot(xi);
            vir += 2*forces[i].dot(xc);
        }
        double beta = betaN*nBeads;
        double En_vir = 1.0/beta + pcP1.getLastEnergy() + 1.0/2.0*vir;
        return En_vir;
    }
}

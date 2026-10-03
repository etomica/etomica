package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPICentVirTIA extends DataSourceScalar {
    protected Box box;
    protected final PotentialComputeField pcP1, pcP1EnTIA;
    protected double betaN;
    protected int nBeads;

    public MeterPICentVirTIA(PotentialComputeField pcP1EnTIA, PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
        super("Stuff", Null.DIMENSION);
        this.pcP1 = pcP1;
        this.pcP1EnTIA = pcP1EnTIA;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;
    }

    @Override
    public double getDataAsScalar() {
        pcP1.computeAll(true);
        pcP1EnTIA.computeAll(false);
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
            vir += forces[i].dot(xc);
        }
        double beta = betaN*nBeads;
        double En_vir = 1.0/2.0/beta + pcP1EnTIA.getLastEnergy() + 1.0/2.0*vir;
//        System.out.println("vir: " +En_vir);
        return En_vir;
    }
}

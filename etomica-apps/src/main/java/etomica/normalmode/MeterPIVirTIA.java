package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIVirTIA extends DataSourceScalar {
    protected Box box;
    protected final PotentialComputeField pcP1, pcP1EnTIA;
    protected double betaN;
    protected int nBeads;

    public MeterPIVirTIA(PotentialComputeField pcP1EnTIA, PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
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
        for (int i = 0; i < nBeads; i++){
            Vector xi = box.getLeafList().get(i).getPosition();
            vir -= forces[i].dot(xi);
        }
        double En_vir = pcP1EnTIA.getLastEnergy() + 1.0/2.0*vir;
//        System.out.println("vir: " +En_vir);
        return En_vir;
    }
}

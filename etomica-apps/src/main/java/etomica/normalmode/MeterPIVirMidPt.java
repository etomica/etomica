package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIVirMidPt extends DataSourceScalar {
    protected Box box;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;

    public MeterPIVirMidPt(PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
        super("Stuff", Null.DIMENSION);
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;
    }

    @Override
    public double getDataAsScalar() {
        Vector[] xbar = new Vector[nBeads];
        Vector[] Fbar = new Vector[nBeads];
        double[] Ubar = new double[nBeads];

        for (int i = 0; i < nBeads; i++){
            xbar[i] = box.getSpace().makeVector();
            Fbar[i] = box.getSpace().makeVector();
            Ubar[i] = 0;
            Vector xi = box.getLeafList().get(i).getPosition();
            Vector xip1;
            if (i == nBeads-1){
                xip1 = box.getLeafList().get(0).getPosition();
            } else {
                xip1 = box.getLeafList().get(i+1).getPosition();
            }
            xbar[i].PE(xi);
            xbar[i].PE(xip1);
            xbar[i].TE(0.5);
            double k2 = 1, k4=24;
            Ubar[i] = 0.5*k2*xbar[i].squared() + 1.0/24.0*k4*xbar[i].squared()*xbar[i].squared();
            Fbar[i].PEa1Tv1(-k2, xbar[i]);
            Fbar[i].PEa1Tv1(-k4/6.0*xbar[i].squared(), xbar[i]);
        }

        double En_vir_bar = 0;
        for (int i = 0; i < nBeads; i++){
            En_vir_bar += Ubar[i] - 0.5*Fbar[i].dot(xbar[i]);
        }
        return En_vir_bar/nBeads;
    }
}

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMAvir extends DataSourceScalar {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;
    protected double[] gk;
    protected double[][] M;
    protected Box box;

    public MeterPIHMAvir(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads, double omega2, Box box, double hbar) {
        super("Stuff", Null.DIMENSION);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;

        gk = new double[nBeads];
        int nK = nBeads/2;
        double beta = this.betaN*nBeads;
        double a2 = (this.betaN*hbar/2.0)*(this.betaN*hbar/2.0)*omega2;
        gk[nK] = -1.0/2.0/beta;
        for(int k = 1; k <= nK; k++){
            double sin = Math.sin(Math.PI*k/nBeads);
            gk[nK-k] = 1.0/2.0/beta*(sin*sin - a2)/(sin*sin + a2);
            if (k != nK || nBeads % 2 != 0){ //odd
                gk[nK+k] = gk[nK-k];
            }
        }

        M = new double[nBeads][nBeads];

        for (int i = 0; i < nBeads ; i++){
            for (int j = 0; j < nBeads ; j++){
                M[i][j] = gk[nK]/nBeads;
                for (int k = 1; k <= nK ; k++){
                    M[i][j] += 1.0/nBeads*gk[nK-k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    if (k != nK || nBeads % 2 != 0){ //odd
                        M[i][j] += 1.0/nBeads*gk[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    }
                }
            }
        }

    }

    @Override
    public double getDataAsScalar() {
        pmBonding.computeAll(true);
        pcP1.computeAll(true);

        double En_HMA =  pcP1.getLastEnergy() ;


        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();
        double beta = betaN*nBeads;

        for (int i = 0; i < nBeads; i++){
            Vector xi = box.getLeafList().get(i).getPosition();
            En_HMA -= (gk[i] + 0.5*forcesU[i].dot(xi));
            for (int j = 0; j < nBeads; j++){
                Vector xj = box.getLeafList().get(j).getPosition();
                if (nBeads == 1) {
                    En_HMA -= beta*(forcesU[i].dot(xj))*M[i][j];
                } else {
                    En_HMA -= beta*(forcesU[i].dot(xj) + forcesK[i].dot(xj))*M[i][j];
                }
            }
        }
//        System.out.println("hma: " + En_HMA);
        return En_HMA;
    }
}

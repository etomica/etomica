package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.util.Constants;

public class MeterPIHMA extends DataSourceScalar {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;
    protected double[] gk;
    protected double[][] M;
    protected Box box;
    public static final double kB = 1.0; // Constants.BOLTZMANN_K;
    public static final double hbar = 1.0;// Constants.PLANCK_H/(2.0*Math.PI);

    public MeterPIHMA(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads, double omega, Box box) {
        super("Stuff", Null.DIMENSION);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;

        gk = new double[nBeads];
        int nK = nBeads/2;
        double beta = this.betaN*nBeads;
        double a = this.betaN*hbar*omega/2.0;
        gk[nK] = -1.0/2.0/beta;
        for(int k = 1; k <= nK; k++){
            double sin = Math.sin(Math.PI*k/nBeads);
            gk[nK-k] = 1.0/2.0/beta*(sin*sin - a*a)/(sin*sin + a*a);
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

        double En_har = 1.0/2.0/this.betaN;
        for (int k = 0; k < nBeads; k++){
            En_har -= gk[k];
        }
        System.out.println(" En_harm: " + En_har);
    }

    @Override
    public double getDataAsScalar() {
        pmBonding.computeAll(true);
        pcP1.computeAll(true);

        double En_HMA = 1.0/2.0/betaN + pcP1.getLastEnergy() - pmBonding.getLastEnergy();


        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();
        double beta = betaN*nBeads;

        for (int i = 0; i < nBeads; i++){
            En_HMA -= gk[i];
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

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMA extends DataSourceScalar {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;
    protected double[] gk;
    protected double[][] M;
    protected Box box;


    public static final double kB = 1.0;//Constants.BOLTZMANN_K


    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);

    public MeterPIHMA(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads, double omega, Box box) {
        super("Stuff", Null.DIMENSION);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;

        gk = new double[nBeads];
        int nK = nBeads/2;
        double beta = betaN*nBeads;
        double a = betaN*hbar*omega/2.0;
        gk[nK] = -1.0/2.0/beta;
        for(int k = 1; k <= nK; k++){
            double sin = Math.sin(Math.PI*k/nBeads);
            gk[nK+k] = 1.0/2.0/beta*(sin*sin - a*a)/(sin*sin + a*a);
            gk[nK-k] = gk[nK+k];
        }

        M = new double[nBeads][nBeads];

        for (int i = 0; i < nBeads ; i++){
            for (int j = 0; j < nBeads ; j++){
                M[i][j] = gk[nK]/nBeads;
                for (int k = 1; k <= nK ; k++){
                    M[i][j] += 2.0/nBeads*gk[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                }
            }
        }
    }

    @Override
    public double getDataAsScalar() {
        pmBonding.computeAll(true);
        pcP1.computeAll(true);

        double En_HMA = 1.0/2.0/betaN + pcP1.getLastEnergy()/nBeads - pmBonding.getLastEnergy();


        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();
        double beta = betaN*nBeads;

        for (int i = 0; i < nBeads; i++){
            En_HMA -= gk[i];
            for (int j = 0; j < nBeads; j++){
                Vector xj = box.getLeafList().get(j).getPosition();
                En_HMA -= beta*(forcesU[i].dot(xj)/nBeads + forcesK[i].dot(xj))*M[i][j];
            }
        }
        return En_HMA;
    }
}

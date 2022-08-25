package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import etomica.util.Constants;

public class MeterPIHMAcent extends DataSourceScalar {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN;
    protected int nBeads;
    protected double[] gk;
    protected double[][] M;
    protected Box box;
    public static final double kB = 1.0; // Constants.BOLTZMANN_K;
    public static final double hbar = 1.0;// Constants.PLANCK_H/(2.0*Math.PI);

    public MeterPIHMAcent(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads, double omega2, Box box) {
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
            double sin2 = Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads);
            gk[nK-k] = 1.0/2.0/beta*(sin2 - a2)/(sin2 + a2);
            if (k != nK || nBeads % 2 != 0){ //odd
                gk[nK+k] = gk[nK-k];
            }
        }

//        M = new double[nBeads-1][nBeads-1];

//        M = new double[][]{{0.484312519224854, 0.000000000000000, 0.007689941556444, 0.007689941556444},
//                {-0.007689941556444, 0.476622577668410, -0.007689941556444, 0},
//                {-0.000000000000000, -0.007689941556444, 0.476622577668410, -0.007689941556444},
//                {0.007689941556444, 0.007689941556444, 0.000000000000000, 0.484312519224854}
//        };


        M = new double[][]  {{  0.496273920820615 ,  0.000000000000000 ,  0.002944568771043 ,  0.005131962412844 ,  0.006580258558806 ,  0.007301426598565 ,  0.007301426598565 ,  0.006580258558806 ,  0.005131962412844 ,  0.002944568771043}
 ,{ -0.002944568771043 ,  0.493329352049571 , -0.002944568771043 , -0.000000000000000 ,  0.002187393641801 ,  0.003635689787763 ,  0.004356857827522 ,  0.004356857827522 ,  0.003635689787763 ,  0.002187393641801}
 ,{ -0.002187393641801 , -0.005131962412844 ,  0.491141958407770 , -0.005131962412844 , -0.002187393641801 ,  0.000000000000000 ,  0.001448296145962 ,  0.002169464185721 ,  0.002169464185721 ,  0.001448296145962}
 ,{ -0.001448296145962 , -0.003635689787763 , -0.006580258558806 ,  0.489693662261809 , -0.006580258558806 , -0.003635689787763 , -0.001448296145962 ,  0.000000000000000 ,  0.000721168039759 ,  0.000721168039759}
 ,{ -0.000721168039759 , -0.002169464185721 , -0.004356857827522 , -0.007301426598565 ,  0.488972494222050 , -0.007301426598565 , -0.004356857827522 , -0.002169464185721 , -0.000721168039759 ,  0.000000000000000}
 ,{ -0.000000000000000 , -0.000721168039759 , -0.002169464185721 , -0.004356857827522 , -0.007301426598565 ,  0.488972494222050 , -0.007301426598565 , -0.004356857827521 , -0.002169464185720 , -0.000721168039759}
 ,{  0.000721168039759 ,  0.000721168039759 ,  0.000000000000000 , -0.001448296145962 , -0.003635689787763 , -0.006580258558806 ,  0.489693662261809 , -0.006580258558806 , -0.003635689787763 , -0.001448296145962}
 ,{  0.001448296145962 ,  0.002169464185721 ,  0.002169464185721 ,  0.001448296145962 ,  0.000000000000000 , -0.002187393641801 , -0.005131962412844 ,  0.491141958407770 , -0.005131962412844 , -0.002187393641801}
 ,{  0.002187393641801 ,  0.003635689787763 ,  0.004356857827521 ,  0.004356857827522 ,  0.003635689787763 ,  0.002187393641801 , -0.000000000000000 , -0.002944568771043 ,  0.493329352049571 , -0.002944568771043}
 ,{  0.002944568771043 ,  0.005131962412844 ,  0.006580258558806 ,  0.007301426598565 ,  0.007301426598565 ,  0.006580258558806 ,  0.005131962412844 ,  0.002944568771043 , -0.000000000000000 ,  0.496273920820615}};




//        for (int i = 0; i < nBeads-1 ; i++){
//            for (int j = 0; j < nBeads-1 ; j++){
//                for (int k = 1; k <= nK ; k++){
//                    M[i][j] += 1.0/nBeads*gk[nK-k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
//                    if (k != nK || nBeads % 2 != 0){ //odd
//                        M[i][j] += 1.0/nBeads*gk[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
//                    }
//                }
//            }
//        }

        double En_har = 1.0/2.0/this.betaN;
        for (int k = 0; k < nBeads; k++){
            En_har -= gk[k];
        }
    }

    @Override
    public double getDataAsScalar() {

        Vector xc = box.getSpace().makeVector();
        for (int i = 0; i < nBeads; i++){
            Vector xi = box.getLeafList().get(i).getPosition();
            xc.PE(xi);
        }
        xc.TE(1.0/nBeads);


        pmBonding.computeAll(true);
        pcP1.computeAll(true);

        double En_HMA = 1.0/2.0/betaN + pcP1.getLastEnergy() - pmBonding.getLastEnergy();


        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();
        double beta = betaN*nBeads;

        int nK = nBeads/2;

        for (int i = 0; i < nBeads; i++){
            En_HMA -= gk[i];
            if (nBeads == 1) {
                En_HMA -= beta*gk[nK] * forcesU[i].dot(xc);
            } else {
                En_HMA -= beta*gk[nK]*(forcesU[i].dot(xc) + forcesK[i].dot(xc));
            }
        }
        for (int i = 0; i < nBeads-1; i++){
            for (int j = 0; j < nBeads-1; j++){
                Vector xj = box.getLeafList().get(j).getPosition();
                Vector yj = box.getSpace().makeVector();
                yj.Ev1Mv2(xj, xc);
                if (nBeads == 1) {
                    En_HMA -= beta*(forcesU[i].dot(yj))*M[j][i];
                    En_HMA += beta*(forcesU[nBeads-1].dot(yj))*M[j][i];
                } else {
                    En_HMA -= beta*(forcesU[i].dot(yj) + forcesK[i].dot(yj))*M[j][i];
                    En_HMA += beta*(forcesU[nBeads-1].dot(yj) + forcesK[nBeads-1].dot(yj))*M[j][i];
                }
            }
        }
        return En_HMA;
    }
}

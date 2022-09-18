package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMA implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double betaN, beta, omegan;
    protected int nBeads;
    protected double[] gk, gk2;
    protected double[][] M, M2;
    protected Box box;
    protected double drdotHdrdot;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rdot, rddot;
    protected double EnShift;

    public static final double kB = 1.0; // Constants.BOLTZMANN_K;
    public static final double hbar = 1.0;// Constants.PLANCK_H/(2.0*Math.PI);

    public MeterPIHMA(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double betaN, int nBeads, double omega2, Box box) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;
        omegan = 1.0/hbar/betaN;
        rdot = new Vector[nBeads];
        rddot = new Vector[nBeads];
        for (int i=0; i<nBeads; i++) {
            rdot[i] = this.box.getSpace().makeVector();
            rddot[i] = this.box.getSpace().makeVector();
        }
        gk = new double[nBeads];
        gk2 = new double[nBeads];
        int nK = nBeads/2;
        beta = this.betaN*nBeads;
        double a2 = (this.betaN*hbar/2.0)*(this.betaN*hbar/2.0)*omega2;
        gk[nK] = -1.0/2.0/beta;
        gk2[nK] = 1.0/2.0/beta/beta;
        for(int k = 1; k <= nK; k++){
            double sin = Math.sin(Math.PI*k/nBeads);
            double num = sin*sin - a2;
            double dnum = -2.0*a2/beta;
            double den = sin*sin + a2;
            double dden = 2.0*a2/beta;
            double den2 = den*den;
            gk[nK-k] = 1.0/2.0/beta*num/den;
            gk2[nK-k] = -1.0/2.0/beta/beta*num/den + 1.0/2.0/beta*(den*dnum-num*dden)/den2;
            if (k != nK || nBeads % 2 != 0){ //odd
                gk[nK+k] = gk[nK-k];
                gk2[nK+k] = gk2[nK-k];
            }
        }

        M = new double[nBeads][nBeads];
        M2 = new double[nBeads][nBeads];
        for (int i = 0; i < nBeads ; i++){
            for (int j = 0; j < nBeads ; j++){
                M[i][j] = gk[nK]/nBeads;
                M2[i][j] = (gk2[nK]+gk[nK]*gk[nK])/nBeads;
                for (int k = 1; k <= nK ; k++){
                    M[i][j] += 1.0/nBeads*gk[nK-k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    M2[i][j] += 1.0/nBeads*(gk2[nK-k] + gk[nK-k]*gk[nK-k])*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    if (k != nK || nBeads % 2 != 0){ //odd
                        M[i][j] += 1.0/nBeads*gk[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                        M2[i][j] += 1.0/nBeads*(gk2[nK+k] + gk[nK+k]*gk[nK+k])*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    }
                }
            }
        }

        double En_harm = 1.0/2.0/this.betaN;
        double Cvn_harm = nBeads/2.0/beta/beta;
        for (int k = 0; k < nBeads; k++){
            En_harm -= gk[k];
            Cvn_harm += gk2[k];
        }
        Cvn_harm *= kB*beta*beta;
        System.out.println(" En_harm: " + En_harm + " Cvn_harm: " + Cvn_harm);
        this.EnShift = 0;
    }

    @Override
    public IData getData() {
        drdotHdrdot = 0 ;
        double[] x = data.getData();
//        System.out.println("******** HMA *************");
        for (int i = 0; i < nBeads; i++) {
            rdot[i].E(0);
            rddot[i].E(0);
//            System.out.println(box.getLeafList().get(i).getPosition());
            for (int j = 0; j < nBeads; j++) {
                Vector rj = box.getLeafList().get(j).getPosition();
                rdot[i].PEa1Tv1(M[i][j], rj);
                rddot[i].PEa1Tv1(M2[i][j], rj);
            }
        }

        pmBonding.computeAll(true);
        pcP1.computeAll(true, this);

        double En = 1.0/2.0/betaN + pcP1.getLastEnergy() - pmBonding.getLastEnergy() - EnShift;
        double Cvn = nBeads/2.0/beta/beta - 2.0*pmBonding.getLastEnergy()/beta;

        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();

        for (int i = 0; i < nBeads; i++){
            En -= gk[i];
            Cvn += gk2[i];
            if (nBeads == 1){
                En -= beta*(forcesU[i].dot(rdot[i]));
                Cvn += 2.0*(forcesU[i].dot(rdot[i]));//rdot
                Cvn += beta*(forcesU[i].dot(rddot[i]));//rddot
            } else {
                En -= beta*(forcesU[i].dot(rdot[i]) + forcesK[i].dot(rdot[i]));
                Cvn += 2.0*(forcesU[i].dot(rdot[i]) - forcesK[i].dot(rdot[i]));//rdot
                Cvn += beta*(forcesU[i].dot(rddot[i]) + forcesK[i].dot(rddot[i]));//rddot
            }
            int ip = i == nBeads-1 ?  0 : i+1;
            Vector tmpV = box.getSpace().makeVector();
            tmpV.Ev1Mv2(rdot[i], rdot[ip]);
            Cvn -= betaN*omegan*omegan*(tmpV.squared());
        }
        Cvn -= beta*drdotHdrdot;
//        System.out.println("hma: " + En);
        x[0] = En;
        x[1] = Cvn;

        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector tmpV = box.getSpace().makeVector();
        tmpV.E(rdot[j]);
        Hij.transform(tmpV);
        drdotHdrdot += rdot[i].dot(tmpV); //drdot.H.drdot
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public boolean wantsHessian() {
        return true;
    }

    public void setShift(double EnShift){
        this.EnShift = EnShift;
    }

}

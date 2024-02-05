package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMA implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1;
    protected double betaN, beta, omega2, omegaN;
    protected int nBeads;
    protected double[] gk, gk2;
    protected double[][] M, M2;
    protected Box box;
    protected double drdotHdrdot;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rdot, rddot;
    protected final Vector[] latticePositions;
    protected int dim, numAtoms;
    protected double EnShift;

    public MeterPIHMA(PotentialMasterBonding pmBonding, PotentialCompute pcP1, double betaN, int nBeads, double omega2, Box box, double hbar) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.dim = box.getSpace().D();
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.box = box;
        this.beta = betaN*nBeads;
        this.omega2 = omega2;
        omegaN = Math.sqrt(nBeads)/(beta*hbar);
        numAtoms = box.getMoleculeList().size();
        rdot = box.getSpace().makeVectorArray(numAtoms*nBeads);
        rddot = box.getSpace().makeVectorArray(numAtoms*nBeads);
        gk = new double[nBeads];
        gk2 = new double[nBeads];
        int nK = nBeads/2;
        beta = betaN*nBeads;
        double a2 = (betaN*hbar/2.0)*(betaN*hbar/2.0)*omega2;
        if (omega2 == 0){
            gk[nK] = 0;
            gk2[nK] = 0;
        } else {
            gk[nK] = -1.0/2.0/beta;
            gk2[nK] = 1.0/2.0/beta/beta;
        }
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

        if (omega2 != 0) {
            double En_ho_nm = dim*nBeads / 2.0 / this.beta;
            double Cvn_ho_nm = dim*nBeads / 2.0 / beta / beta;
            for (int k = 0; k < nBeads; k++) {
                En_ho_nm -= dim*gk[k];
                Cvn_ho_nm += dim*gk2[k];
            }
            Cvn_ho_nm *= beta * beta;
            //COM
            En_ho_nm -= dim/2.0/beta;
            Cvn_ho_nm -= dim/2.0/beta/beta;
            System.out.println(" En_ho_nm:  " + En_ho_nm);
            System.out.println(" Cvn_ho_nm: " + Cvn_ho_nm);
        }


        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }

        this.EnShift = 0;
    }

    @Override
    public IData getData() {
        drdotHdrdot = 0 ;
        double[] x = data.getData();
        Vector drShift = box.getSpace().makeVector();
        if (box.getMoleculeList().size() > 1) {
            drShift = computeShift();
        }
        IMoleculeList molecules = box.getMoleculeList();
        Vector drj = box.getSpace().makeVector();
        for (IMolecule m : molecules) {
            IAtomList atoms = m.getChildList();
            for (int i = 0; i < nBeads; i++) {
                int ia = atoms.get(i).getLeafIndex();
                rdot[ia].E(0);
                rddot[ia].E(0);
            }

            for (int j = 0; j < nBeads; j++) {
                computeDR(atoms.get(j).getPosition(), latticePositions[m.getIndex()], drShift, drj);
                for (int i=0; i<nBeads; i++) {
                    int ia = atoms.get(i).getLeafIndex();
                    rdot[ia].PEa1Tv1(M[i][j], drj);
                    rddot[ia].PEa1Tv1(M2[i][j], drj);
                }
            }
        }

        pmBonding.computeAll(true);
//        pcP1.computeAll(true, null); //no Cv
        pcP1.computeAll(true, this); //with Cv (replace 'this' by 'null' for no Cv)

        double En = dim*numAtoms*nBeads/2.0/beta + pcP1.getLastEnergy() - pmBonding.getLastEnergy();
        double Cvn = dim*numAtoms*nBeads/2.0/beta/beta - 2.0/beta*pmBonding.getLastEnergy();
        for (int i=0; i<nBeads; i++) {
            En -= dim*numAtoms*gk[i];
            Cvn += dim*numAtoms*gk2[i];
        }
        //COM
        if (numAtoms > 1 && omega2 != 0) {
            En -= dim/2.0/beta; //com: -dim*gk[nK]
            Cvn -= dim/2.0/beta/beta;
        }


        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();

        double massRing = box.getLeafList().get(0).getType().getMass()*nBeads;
        Vector tmpV = box.getSpace().makeVector();

        for (IMolecule molecule : box.getMoleculeList()) {
            IAtomList atoms = molecule.getChildList();
            for (int jj=0; jj<atoms.size(); jj++) {
                int j = atoms.get(jj).getLeafIndex();
                En -= beta*(forcesU[j].dot(rdot[j]) + forcesK[j].dot(rdot[j]));
                Cvn += 2.0*(forcesU[j].dot(rdot[j]) - forcesK[j].dot(rdot[j]));
                Cvn += beta*(forcesU[j].dot(rddot[j]) + forcesK[j].dot(rddot[j]));
                int jjp = (jj+1)%nBeads;
                int jp = atoms.get(jjp).getLeafIndex();
                tmpV.Ev1Mv2(rdot[j], rdot[jp]);
                Cvn -= beta*massRing*omegaN*omegaN*(tmpV.squared());
            }
        }

        Cvn -= beta*drdotHdrdot;
        x[0] = En - EnShift;
        x[1] = Cvn + (En - EnShift)*(En - EnShift);

        return data;
    }

    public void fieldComputeHessian(int i, Tensor Hii) {
        Vector tmpV = box.getSpace().makeVector();
        tmpV.E(rdot[i]);
        Hii.transform(tmpV);
        drdotHdrdot += rdot[i].dot(tmpV);
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector rdotij = box.getSpace().makeVector();
        rdotij.Ev1Mv2(rdot[i], rdot[j]);
        Hij.transform(rdotij);
        drdotHdrdot += rdotij.dot(rdot[i]) - rdotij.dot(rdot[j]);
    }

    protected Vector computeShift() {
        if (box.getMoleculeList().size() == 1) {
            return box.getSpace().makeVector();
        }
        int n = box.getMoleculeList().get(0).getChildList().size();
        Vector shift0 = box.getSpace().makeVector();
        Boundary boundary = box.getBoundary();
        Vector dr = box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        dr.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[0]);
        boundary.nearestImage(dr);
        shift0.Ea1Tv1(-1, dr);
        // will shift ring0 back to lattice site; everything should be close and PBC should lock in
        // now determine additional shift needed to bring back to original COM
        Vector totalShift = box.getSpace().makeVector();
        for (int j = 0; j < box.getMoleculeList().size(); j++) {
            IMolecule m = box.getMoleculeList().get(j);
            for (int i = 0; i < n; i++) {
                Vector r = m.getChildList().get(i).getPosition();
                dr.Ev1Mv2(r, latticePositions[j]);
                dr.PE(shift0);
                boundary.nearestImage(dr);
                totalShift.PE(dr);
            }
        }
        totalShift.TE(-1.0/box.getLeafList().size());
        totalShift.PE(shift0);
        return totalShift;
    }

    protected void computeDR(Vector r, Vector latticeSite, Vector shift, Vector dr) {
        dr.Ev1Mv2(r, latticeSite);
        dr.PE(shift);
        box.getBoundary().nearestImage(dr);
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

    public void setEnShift(double E) { this.EnShift = E; }
}

package etomica.normalmode;

import etomica.atom.IAtom;
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

public class MeterPIHMAInf implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1harm, pcP1ah;
    protected double temperature, beta, omegaN, R, sinhR, cothR;
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
    protected int dim;
    protected double fac1, fac2, fac3, fac4, mOmega1, mOmega2;

    public MeterPIHMAInf(PotentialMasterBonding pmBonding, PotentialCompute pcP1harm,PotentialCompute pcP1ah, double temperature, int nBeads, double omega,double omegaSample, Box box, double hbar) {
        int nData = 2;
        this.box = box;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.dim = box.getSpace().D();
        this.pmBonding = pmBonding;
        this.pcP1harm = pcP1harm;
        this.pcP1ah = pcP1ah;
        this.temperature = temperature;
        this.nBeads = nBeads;
        this.beta = 1/temperature;
        this.omegaN = Math.sqrt(nBeads)/(beta*hbar);
        double mass = box.getLeafList().get(0).getType().getMass()*nBeads;

        double RSample = beta*hbar*Math.sqrt(omegaSample)/nBeads;
        int nAtoms = box.getLeafList().size();
        rdot = box.getSpace().makeVectorArray(nAtoms);
        rddot = box.getSpace().makeVectorArray(nAtoms);
        gk = new double[nBeads];
        gk2 = new double[nBeads];
        int nK = nBeads/2;
        if (omegaSample == 0){
            gk[nK] = 0;
            gk2[nK] = 0;
        } else {
            gk[nK] = -RSample/2.0/beta/Math.sinh(RSample);
            gk2[nK] = RSample*RSample/2/beta/beta*Math.cosh(RSample)/Math.sinh(RSample)/Math.sinh(RSample);
        }
        double dRSample = RSample/beta;
        for(int k = 1; k <= nK; k++){
            double sin = Math.sin(Math.PI*k/nBeads);
            double num = 2/Math.tanh(RSample)*sin*sin - Math.tanh(RSample/2);
            double dnum = -2/Math.sinh(RSample)/Math.sinh(RSample)*sin*sin-1.0/2.0/Math.cosh(RSample/2)/Math.cosh(RSample/2);
            double den = sin*sin + Math.sinh(RSample/2)*Math.sinh(RSample/2);
            double dden = Math.sinh(RSample)/2;
            double den2 = den*den;
            if (omegaSample == 0) {
                gk[nK-k] = 1.0/2.0/beta;
                gk2[nK-k] = -1.0/2.0/beta/beta;
            } else {
                gk[nK-k] =  RSample/4.0/beta*num/den;
                gk2[nK-k] = RSample/4/beta*(den*dnum-num*dden)/den2*dRSample;
            }

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

        this.R = beta*hbar*omega/nBeads;
        sinhR = Math.sinh(R);
        cothR = 1/Math.tanh(R);
        this.fac1 = R*R/beta*(1-2*cothR*cothR);
        this.fac2 = R*R/2/beta/Math.cosh(R/2)/Math.cosh(R/2);
        this.fac3 = -R*cothR;
        this.fac4 = R/sinhR;
        this.mOmega1 = mass*omega*omega/nBeads/R/sinhR;
        this.mOmega2 = 2*mass*omega*omega*Math.tanh(R/2)/nBeads/R;

        if (omegaSample != 0) {
            double En_ho_nm = hbar*omega/2.0*cothR;
            double Cvn_ho_nm = nBeads / 2.0 / beta / beta;
            for (int k = 0; k < nBeads; k++) {
                En_ho_nm -= gk[k];
                Cvn_ho_nm += gk2[k];
            }
            Cvn_ho_nm *= beta * beta;
            System.out.println(" En_ho_nm:  " + En_ho_nm);
            System.out.println(" Cvn_ho_nm: " + Cvn_ho_nm);
        }


        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }
    }

    @Override
    public IData getData() {
        drdotHdrdot = 0 ;
        double[] x = data.getData();
        Vector shift = computeShift();
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
                computeDR(atoms.get(j).getPosition(), latticePositions[m.getIndex()], shift, drj);
                for (int i=0; i<nBeads; i++) {
                    int ia = atoms.get(i).getLeafIndex();
                    rdot[ia].PEa1Tv1(M[i][j], drj);
                    rddot[ia].PEa1Tv1(M2[i][j], drj);
                }
            }
        }

        pmBonding.computeAll(true);
        pcP1harm.computeAll(true); //no Cv
//        pcP1ah.computeAll(true); //no Cv
        pcP1ah.computeAll(true, this); //with Cv (replace 'this' by 'null' for no Cv)

        int numAtoms = molecules.size();
        double En = dim*nBeads*numAtoms/2.0/beta*R*cothR + R/sinhR*pcP1harm.getLastEnergy() + pcP1ah.getLastEnergy() - R*cothR*pmBonding.getLastEnergy();
        double Cvn = dim*nBeads*numAtoms/2.0/beta/beta*R*R/sinhR/sinhR + fac1*pmBonding.getLastEnergy() + fac2*pcP1harm.getLastEnergy();

        if (box.getMoleculeList().size() > 1) {
            En -= dim/2.0/beta;
        }
        for (int i=0; i<nBeads; i++) {
            En -= dim*numAtoms*gk[i];
            Cvn += dim*numAtoms*gk2[i];
        }

        Vector[] forcesUharm = pcP1harm.getForces();
        Vector[] forcesUah = pcP1ah.getForces();
        Vector[] forcesK = pmBonding.getForces();

        for (IAtom a : box.getLeafList()) {
            int j = a.getLeafIndex();

            En -= beta*( forcesK[j].dot(rdot[j])+forcesUharm[j].dot(rdot[j]) + forcesUah[j].dot(rdot[j]));
            Cvn += beta*(forcesK[j].dot(rddot[j])+forcesUharm[j].dot(rddot[j]) + forcesUah[j].dot(rddot[j]));//rddot
            Cvn += 2.0*(fac3*forcesK[j].dot(rdot[j])+fac4*forcesUharm[j].dot(rdot[j]) + forcesUah[j].dot(rdot[j]));//rdot

            int jp = a.getIndex() == nBeads-1 ?  (j-nBeads+1) : j+1;
            Vector tmpV = box.getSpace().makeVector();
            tmpV.Ev1Mv2(rdot[j], rdot[jp]);
            Cvn -= (beta*mOmega1*(tmpV.squared()) + beta*mOmega2*rdot[j].squared());
        }
        Cvn -= beta*drdotHdrdot;
        x[0] = En;
        x[1] = Cvn + x[0]*x[0];

        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector tmpV = box.getSpace().makeVector();
        tmpV.E(rdot[j]);
        Hij.transform(tmpV);
        drdotHdrdot += rdot[i].dot(tmpV); //drdot.H.drdot
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

}

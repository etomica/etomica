package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.molecule.CenterOfMass;
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
    protected final Vector[] latticePositions;

    public MeterPIHMA(PotentialMasterBonding pmBonding, PotentialCompute pcP1, double betaN, int nBeads, double omega2, Box box, double hbar) {
        int nData = 1;
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
        int nAtoms = box.getLeafList().size();
        rdot = box.getSpace().makeVectorArray(nAtoms);
        rddot = box.getSpace().makeVectorArray(nAtoms);
        gk = new double[nBeads];
        gk2 = new double[nBeads];
        int nK = nBeads/2;
        beta = betaN*nBeads;
        double a2 = (betaN*hbar/2.0)*(betaN*hbar/2.0)*omega2;
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
        Cvn_harm *= beta*beta;

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
        Vector shift = computeShift();
//        System.out.println("******** HMA *************");
        IMoleculeList molecules = box.getMoleculeList();
        Vector dr = box.getSpace().makeVector();
        for (IMolecule m : molecules) {
            IAtomList atoms = m.getChildList();
            Vector com = CenterOfMass.position(box, m);
            for (int i = 0; i < nBeads; i++) {
                int ia = atoms.get(i).getLeafIndex();
                rdot[ia].E(0);
                rddot[ia].E(0);
//            System.out.println(box.getLeafList().get(i).getPosition());
                for (int j = 0; j < nBeads; j++) {
                    computeDR(atoms.get(j).getPosition(), com, shift, dr);
                    rdot[ia].PEa1Tv1(M[i][j], dr);
                    rddot[ia].PEa1Tv1(M2[i][j], dr);
                }
            }
        }

        pmBonding.computeAll(true);
        pcP1.computeAll(true, this);

        double En = 1.0/2.0/betaN + pcP1.getLastEnergy() - pmBonding.getLastEnergy() - EnShift;
        double Cvn = nBeads/2.0/beta/beta - 2.0*pmBonding.getLastEnergy()/beta;

        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();

        for (IAtom a : box.getLeafList()) {
            int i = a.getIndex();
            int j = a.getLeafIndex();
            En -= gk[i];
            Cvn += gk2[i];
            if (nBeads == 1){
                En -= beta*(forcesU[j].dot(rdot[j]));
                Cvn += 2.0*(forcesU[j].dot(rdot[j]));//rdot
                Cvn += beta*(forcesU[j].dot(rddot[j]));//rddot
            } else {
                En -= beta*(forcesU[j].dot(rdot[j]) + forcesK[j].dot(rdot[j]));
                Cvn += 2.0*(forcesU[j].dot(rdot[j]) - forcesK[j].dot(rdot[j]));//rdot
                Cvn += beta*(forcesU[j].dot(rddot[j]) + forcesK[j].dot(rddot[j]));//rddot
            }
            int jp = i == nBeads-1 ?  (j-nBeads+1) : j+1;
            Vector tmpV = box.getSpace().makeVector();
            tmpV.Ev1Mv2(rdot[j], rdot[jp]);
            Cvn -= betaN*omegan*omegan*(tmpV.squared());
        }
        Cvn -= beta*drdotHdrdot;
//        System.out.println("hma: " + En);
        x[0] = En;
//        x[1] = Cvn;

        return data;
    }

    protected Vector computeShift() {
        int n = box.getMoleculeList().get(0).getChildList().size();
        Vector shift0 = box.getSpace().makeVector();
        Boundary boundary = box.getBoundary();
        Vector dr = box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        for (int i=0; i<n; i++) {
            dr.Ev1Mv2(atoms.get(i).getPosition(), latticePositions[0]);
            boundary.nearestImage(dr);
            shift0.PE(dr);
        }
        shift0.TE(-1.0/n);
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
        totalShift.TE(-1.0/(n*(box.getMoleculeList().size()-1)));
        totalShift.PE(shift0);
        return totalShift;
    }

    protected void computeDR(Vector r, Vector com, Vector shift, Vector dr) {
        dr.Ev1Mv2(r, com);
        dr.PE(shift);
        box.getBoundary().nearestImage(dr);
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

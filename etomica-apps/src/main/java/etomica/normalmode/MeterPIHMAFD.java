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
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMAFD implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1;
    protected double beta, omega2;
    protected double  omegaN,  omegaN_p,  omegaN_m;
    protected int nBeads;
    protected double[] gk, gk_p, gk_m, gk2;
    protected double[] lambda, lambda_p, lambda_m;
    protected double[][] M, M_p, M_m;
    protected Box box;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rdot;
    protected final Vector[] latticePositions;
    protected int dim, numAtoms;
    protected double EnShift, dbeta, beta_p, beta_m;
    protected double En, dbUdb, hbar;
    protected Vector drShift;
    double[][] eigenvectors;
    protected Vector[] rOrig;

    public MeterPIHMAFD(PotentialMasterBonding pmBonding, PotentialCompute pcP1, double beta, int nBeads, double omega2, Box box, double hbar, double dbeta) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.dim = box.getSpace().D();
        this.pcP1 = pcP1;
        this.pmBonding = pmBonding;
        this.nBeads = nBeads;
        this.box = box;
        this.beta = beta;
        this.dbeta = dbeta;
        beta_p = beta + dbeta;
        beta_m = beta - dbeta;
        numAtoms = box.getMoleculeList().size();
        this.omega2 = omega2;
        this.hbar = hbar;
        rOrig = box.getSpace().makeVectorArray(box.getLeafList().size());
        rdot = box.getSpace().makeVectorArray(numAtoms*nBeads);

        getMappingParams();

        if (omega2 != 0) {
            double En_ho_nm = dim*nBeads/2.0/beta;
            for (int k = 0; k < nBeads; k++) {
                En_ho_nm -= dim*gk[k];
            }
            //COM
            En_ho_nm -= dim/2.0/beta;
            System.out.println(" En_ho_nm:  " + En_ho_nm);
        }

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule molecule : box.getMoleculeList()) {
            latticePositions[molecule.getIndex()].E(molecule.getChildList().get(0).getPosition());
        }

        this.EnShift = 0;
        this.drShift = box.getSpace().makeVector();
    }

    @Override
    public IData getData() {
        double[] x = data.getData();

        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                rOrig[atom.getLeafIndex()].E(atom.getPosition());
            }
        }

        getEn(beta, lambda, M);
        x[0] = En - EnShift;
        if (numAtoms > 1 && omega2 != 0) {
            x[0] -= dim/2.0/beta; //com: -dim*gk[nK]
        }

        double Cvn = dim*numAtoms*nBeads/2.0/beta/beta;
        for (int k = 0; k < nBeads; k++) {
            Cvn += dim*numAtoms*gk2[k];
        }
        if (numAtoms > 1 && omega2 != 0) Cvn -= dim/2.0/beta/beta; //com

        getEn(beta_p, lambda_p, M_p);
        double dbUdb_p = dbUdb;
        getEn(beta_m, lambda_m, M_m);
        double dbUdb_m = dbUdb;
        double d2bUdb2 = (dbUdb_p - dbUdb_m)/(2*dbeta);
//        System.out.println(" FD: " + (-d2bUdb2));
//        System.out.println();
        x[1] = Cvn - d2bUdb2  + x[0]*x[0];

        return data;
    }

    private void getEn(double beta_new, double[] lambda_new, double[][] M_new) {
        if (beta_new != beta) scaleCoord(beta_new, lambda_new);
        if (numAtoms > 1 && beta_new == beta) drShift = computeShift();

        // rdot
        Vector drj = box.getSpace().makeVector();
        for (int i = 0; i < rdot.length; i++) rdot[i].E(0);
        for (IMolecule molecule : box.getMoleculeList()) {
            IAtomList atoms = molecule.getChildList();
            for (int j = 0; j < nBeads; j++) {
                drj.Ev1Mv2(atoms.get(j).getPosition(), latticePositions[molecule.getIndex()]);
                drj.PE(drShift);
                box.getBoundary().nearestImage(drj);
                for (int i = 0; i < nBeads; i++) {
                    int iAtom = atoms.get(i).getLeafIndex();
                    rdot[iAtom].PEa1Tv1(M_new[i][j], drj);
                }
            }
        }

        pcP1.computeAll(true);
        pmBonding.computeAll(true);
        double beta2 = beta*beta;
        double beta_new2 = beta_new*beta_new;
        double U = pcP1.getLastEnergy();
        double K = beta2/beta_new2*pmBonding.getLastEnergy();
        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();
        for (int i = 0; i < forcesK.length; i++) {
            forcesK[i].TE(beta2/beta_new2);
        }

        dbUdb = U - K;
        En = dim*numAtoms*nBeads/2.0/beta + U - K;
        for (int i = 0; i < nBeads; i++) En -= dim*numAtoms*gk[i];

        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                int j = atom.getLeafIndex();
                dbUdb -= beta_new*(forcesU[j].dot(rdot[j]) + forcesK[j].dot(rdot[j]));
                En -= beta_new*(forcesU[j].dot(rdot[j]) + forcesK[j].dot(rdot[j]));
            }
        }

        // If atoms scaled, bring them back to their original coords
        if (beta_new != beta) {
            for (IMolecule molecule : box.getMoleculeList()) {
                for (IAtom atom : molecule.getChildList()) {
                    atom.getPosition().E(rOrig[atom.getLeafIndex()]);
                }
            }
        }

    }

    // Scale coordinates using NM
    private void scaleCoord(double beta_new, double[] lambda_new) {
        if (numAtoms > 1) drShift = computeShift();
        Vector dri = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            // Real -> NM
            Vector[] q = box.getSpace().makeVectorArray(nBeads);
            for (IAtom atom : molecule.getChildList()) {
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                box.getBoundary().nearestImage(dri);
                for (int k = 0; k < nBeads; k++) {
                    q[k].PEa1Tv1(eigenvectors[atom.getIndex()][k]/Math.sqrt(nBeads), dri);
                }
            }
            // Scale NM
            for (int k = 0; k < nBeads; k++) {
                q[k].TE(Math.sqrt(lambda[k]*beta/lambda_new[k]/beta_new));
            }
            // Scaled NM -> Real
            for (IAtom atom : molecule.getChildList()) {
                dri = box.getSpace().makeVector();
                for (int k = 0; k < nBeads; k++) {
                    dri.PEa1Tv1(Math.sqrt(nBeads)*eigenvectors[atom.getIndex()][k], q[k]);
                }
                box.getBoundary().nearestImage(dri);
                atom.getPosition().Ev1Pv2(latticePositions[molecule.getIndex()], dri);
            }
        }
    }

    // gk, gk2, eigenvalues and M at 0,+,- perturpations.
    // eigenvectors
    private void getMappingParams() {
        int nK = nBeads/2;
        //Scaling parameters: gk
        gk = new double[nBeads]; gk_p = new double[nBeads]; gk_m = new double[nBeads];
        gk2 = new double[nBeads];
        if (omega2 == 0){
            gk[nK] = 0; gk_p[nK] = 0; gk_m[nK] = 0;
            gk2[nK] = 0;
        } else {
            gk[nK] = -1.0/2.0/beta; gk_p[nK] = -1.0/2.0/beta_p; gk_m[nK] = -1.0/2.0/beta_m;
            gk2[nK] = 1.0/2.0/beta/beta;
        }
        double a2 = Math.pow(beta*hbar/nBeads/2.0,2.0)*omega2;
        double a2_p = Math.pow(beta_p*hbar/nBeads/2.0,2.0)*omega2;
        double a2_m = Math.pow(beta_m*hbar/nBeads/2.0,2.0)*omega2;
        for(int k = 1; k <= nK; k++){
            double sin = Math.sin(Math.PI*k/nBeads);
            double num = sin*sin - a2 , den = sin*sin + a2;
            double num_p = sin*sin - a2_p , den_p = sin*sin + a2_p;
            double num_m = sin*sin - a2_m , den_m = sin*sin + a2_m;
            double dnum = -2.0*a2/beta, dden = 2.0*a2/beta;
            double den2 = den*den;
            gk[nK-k] = 1.0/2.0/beta*num/den;  gk_p[nK-k] = 1.0/2.0/beta_p*num_p/den_p;  gk_m[nK-k] = 1.0/2.0/beta_m*num_m/den_m;
            gk2[nK-k] = -1.0/2.0/beta/beta*num/den + 1.0/2.0/beta*(den*dnum-num*dden)/den2;
            if (k != nK || nBeads % 2 != 0){ //odd
                gk[nK+k] = gk[nK-k];  gk_p[nK+k] = gk_p[nK-k];  gk_m[nK+k] = gk_m[nK-k];
                gk2[nK+k] = gk2[nK-k];
            }
        }

        //Mapping matrix
        M = new double[nBeads][nBeads];  M_p = new double[nBeads][nBeads];  M_m = new double[nBeads][nBeads];
        for (int i = 0; i < nBeads ; i++){
            for (int j = 0; j < nBeads ; j++){
                M[i][j] = gk[nK]/nBeads;  M_p[i][j] = gk_p[nK]/nBeads;  M_m[i][j] = gk_m[nK]/nBeads;
                for (int k = 1; k <= nK ; k++){
                    M[i][j] += 1.0/nBeads*gk[nK-k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    M_p[i][j] += 1.0/nBeads*gk_p[nK-k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    M_m[i][j] += 1.0/nBeads*gk_m[nK-k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    if (k != nK || nBeads % 2 != 0){ //odd
                        M[i][j] += 1.0/nBeads*gk[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                        M_p[i][j] += 1.0/nBeads*gk_p[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                        M_m[i][j] += 1.0/nBeads*gk_m[nK+k]*Math.cos(2.0*Math.PI*k/nBeads*(i-j));
                    }
                }
            }
        }

        // Eigenvalues
        double mass = box.getLeafList().get(0).getType().getMass();
        lambda = new double[nBeads];  lambda_p = new double[nBeads];  lambda_m = new double[nBeads];
        lambda[nK] = mass*omega2; lambda_p[nK] = mass*omega2; lambda_m[nK] = mass*omega2;
        omegaN = nBeads/(beta*hbar);  omegaN_p = nBeads/(beta_p*hbar);  omegaN_m = nBeads/(beta_m*hbar);
        double lambda_k;
        for(int k = 1; k <= (nBeads-1)/2; k++){
            lambda_k = 4.0*mass*omegaN*omegaN*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda[nK-k] = lambda_k; lambda[nK+k] = lambda_k;
            lambda_k = 4.0*mass*omegaN_p*omegaN_p*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda_p[nK-k] = lambda_k; lambda_p[nK+k] = lambda_k;
            lambda_k = 4.0*mass*omegaN_m*omegaN_m*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda_m[nK-k] = lambda_k; lambda_m[nK+k] = lambda_k;
        }
        if (nBeads % 2 == 0){
            int k = nK;
            lambda_k = 4.0*mass*omegaN*omegaN*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda[0] = lambda_k;
            lambda_k = 4.0*mass*omegaN_p*omegaN_p*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda_p[0] = lambda_k;
            lambda_k = 4.0*mass*omegaN_m*omegaN_m*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda_m[0] = lambda_k;
        }

        // Eigenvectors
        eigenvectors = new double[nBeads][nBeads];
        for (int i = 0; i < nBeads; i++) {
            eigenvectors[i][nK] = 1.0/Math.sqrt(nBeads);//k=0
            for (int k = 1; k <= (nBeads-1)/2; k++) {
                double arg = 2.0*Math.PI/nBeads*i*k;
                eigenvectors[i][nK-k] = Math.sqrt(2.0)*Math.sin(-arg)/Math.sqrt(nBeads);
                eigenvectors[i][nK+k] = Math.sqrt(2.0)*Math.cos(arg)/Math.sqrt(nBeads);
            }
            if (nBeads % 2 == 0){ //even
                eigenvectors[i][0] =  Math.pow(-1, i)/Math.sqrt(nBeads);
            }
        }
    }

    protected Vector computeShift() {
        if (numAtoms == 1) return box.getSpace().makeVector();
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
        for (int j = 0; j < numAtoms; j++) {
            IMolecule molecule = box.getMoleculeList().get(j);
            for (int i = 0; i < n; i++) {
                Vector r = molecule.getChildList().get(i).getPosition();
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

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setEnShift(double E) { this.EnShift = E; }
}

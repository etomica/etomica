package etomica.normalmode;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;


public class MeterFM implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pc;
    protected final PotentialCompute pcSS;
    protected double beta;
    protected double trH;
    protected double sumrHr;
    protected double[][] arrayH;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected int numAtoms;
    protected int dof;
    protected int step = 0;
    protected double alpha, dbeta;
    protected Vector[] pos0, lattice, xmin, rmin;
    protected Vector[] rOld;
    protected double avgU, avgF2, avgH, avg4;
    protected double Ulat;
    protected  Vector[] Fold;

    public MeterFM(PotentialCompute pc, PotentialCompute pcSS, double temperature, Box box, double alpha, double dbeta) {
        int nData = 1;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("FM",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pc = pc;
        this.pcSS = pcSS;
        this.beta = 1/temperature;
        this.box = box;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        dof = dim*numAtoms;
        pos0 = box.getSpace().makeVectorArray(numAtoms);
        rOld = box.getSpace().makeVectorArray(numAtoms);
        lattice = box.getSpace().makeVectorArray(numAtoms);
        xmin = box.getSpace().makeVectorArray(numAtoms);
        rmin = box.getSpace().makeVectorArray(numAtoms);

        for (IAtom atom : box.getLeafList()) {
            lattice[atom.getLeafIndex()].E(atom.getPosition());
        }

        this.avgU = 0;
        this.avgF2 = 0;
        this.avgH = 0;
        this.avg4 = 0;
        this.Ulat = 0;
        this.alpha = alpha;
        this.dbeta = dbeta;

        Fold = box.getSpace().makeVectorArray(numAtoms);
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        IAtomList atoms = box.getLeafList();
        for (IAtom atom : atoms) {
            pos0[atom.getLeafIndex()].E(atom.getPosition());
        }

        double beta1 = beta + dbeta;
        arrayH = new double[dof][dof];
        double[][] arrayF = new double[dof][dof];

        trH = 0;
        sumrHr = 0;
        pc.computeAll(false);
//        pc.computeAll(true, this);
        double U0 = pc.getLastEnergy();
        x[0] = U0/numAtoms;
//        Matrix matrix = new Matrix(arrayH);
//        EigenvalueDecomposition ed = matrix.eig();
//        double[] eVals = ed.getRealEigenvalues();
//
//        for (int i=0; i<dof; i++) {
//            if (eVals[i] > 0) System.out.println(Math.sqrt(eVals[i]) + "  p");
//            if (eVals[i] < 0) System.out.println(Math.sqrt(-eVals[i]) + "  m");
//        }

//        double uavg = -4.091361494569530e+00;
//        double usd = 4.285958737466034e-02;
//        double havg = 8.369057047490744e+02;
//        double hsd = 2.303957853643375e+01;
//        System.out.println(step + "  " + (U0/numAtoms-uavg)/usd + "   " + (trH/numAtoms-havg)/hsd);

//        Vector[] F0 = box.getSpace().makeVectorArray(numAtoms);
//        double sumF2 = 0;
//        double sumH = 0;
//        double sumF2H = 0;
//        for (int i = 0; i < numAtoms; i++) {
//            F0[i].E(pc.getForces()[i]);
//            sumF2 += F0[i].squared();
//        }
//        double sumF4 = sumF2*sumF2;
//
//        double sumFHF = 0;
//        for (int i = 0; i < numAtoms; i++){
//            for (int j = 0; j < numAtoms; j++){
//                for (int a = 0; a < dim; a++) {
//                    for (int b = 0; b < dim; b++) {
//                        sumFHF += F0[i].getX(a)*arrayH[dim*i+a][dim*j+b]*F0[j].getX(b);
//                    }
//                }
//            }
//        }
//
//
//        double[][] arrayFF = new double[dof][dof];
//        double[][] arrayFFH = new double[dof][dof];
//
//        for (int i = 0; i < numAtoms; i++) {
//            for (int j = 0; j < numAtoms; j++) {
//                for (int a = 0; a < dim; a++) {
//                    for (int b = 0; b < dim; b++) {
//                        arrayFF[dim*i+a][dim*j+b] = F0[i].getX(a)*F0[j].getX(b);
//                    }
//                }
//            }
//        }
//
//        double sumH2 = 0;
//        for (int i = 0; i < dof; i++) {
//            sumH += arrayH[i][i];
//            for (int j = 0; j < dof; j++) {
//                sumH2 += arrayH[i][j]*arrayH[i][j];
//                for (int k = 0; k < dof; k++) {
//                    arrayFFH[i][j] += arrayFF[i][k]*arrayH[k][j];
//                }
//            }
//        }
//
//        Matrix matrixH = new Matrix(arrayH);
//        Matrix matrixFF = new Matrix(arrayFF);
//        Matrix matrixFFH = new Matrix(arrayFFH);
//
////        double u0avg = -4.094024807648186;//avg @ r=0.8, T=1
//        double u0avg = -4;
//        double dU0 = U0 - u0avg*numAtoms;
//        double scale = 1.0;
//        double c1 = 2*scale*dbeta/beta1*dU0/sumF4;
//        double c2 = -scale*dbeta/beta1/sumF2;
//        double c3 = -scale*dbeta/beta1*dU0/sumF2;
//        matrixFFH = matrixFFH.times(c1);
//        matrixFF = matrixFF.times(c2);
////        matrixH = matrixH.times(c3);
//
//        alpha = 1.448222044961031e-04;
//        matrixH = matrixH.times(-alpha);
//
//        Matrix matrixJ = Matrix.identity(dof,dof);
////        matrixJ.plusEquals(matrixFFH);
////        matrixJ.plusEquals(matrixFF);
//        matrixJ.plusEquals(matrixH);
//        double detJ = matrixJ.det();
//
////        alpha = scale*dbeta*dU0/beta1/sumF2;
//        for (IAtom atom : box.getLeafList()) {
//            atom.getPosition().PEa1Tv1(alpha, F0[atom.getLeafIndex()]);
//            Vector shift = box.getBoundary().centralImage(atom.getPosition());
//            atom.getPosition().PE(shift);
//        }
//        pc.init();
//        pc.computeAll(false);
//        double U = pc.getLastEnergy();
//        double bU0 = beta*U0;
//        double bU = beta1*U;
//        double phi0 = dbeta*U0/numAtoms;
//        double phi = (bU - bU0 - Math.log(detJ))/numAtoms;
//        double phia = (dbeta*U0 - beta1*alpha*sumF2 + alpha*sumH)/numAtoms;
//        double phia2 = (dbeta*U0 - beta1*alpha*sumF2 + alpha*sumH + beta1*0.5*alpha*alpha*sumFHF + 0.5*alpha*alpha*sumH2)/numAtoms;
//
////        double y = (beta1*sumF2 - sumH)/(beta1*sumFHF+sumH2);
////        System.out.println(y);
//        System.out.println(step + "  " + phi0 + "  " + phi + "  " + (phi-phi0));
////        System.out.println(step + "  " + phi0 + "  " + phi + "  " + phia + " " + phia2);
//
//        for (IAtom atom : box.getLeafList()) {
//            atom.getPosition().E(pos0[atom.getLeafIndex()]);
//        }
//        pc.init();
//        pc.computeAll(false);

        step++;
        return data;
    }

    public void fieldComputeHessian(int i, Tensor Hii) { }

    public void pairComputeHessian(int i, int j, Tensor Hij) {
        trH += 2*Hij.trace();
        for (int a = 0; a < dim; a++) {
            for (int b = 0; b < dim; b++) {
                arrayH[dim*i+a][dim*j+b] = -Hij.component(a,b);
                arrayH[dim*j+b][dim*i+a] = -Hij.component(a,b);
                //self
                arrayH[dim*i+a][dim*i+b] += Hij.component(a,b);
                arrayH[dim*j+a][dim*j+b] += Hij.component(a,b);
            }
        }

        Vector rij = box.getSpace().makeVector();
        rij.Ev1Mv2(box.getLeafList().get(i).getPosition() , box.getLeafList().get(j).getPosition());
        box.getBoundary().nearestImage(rij);
        for (int a = 0; a < dim; a++) {
            for (int b = 0; b < dim; b++) {
                sumrHr += rij.getX(a)*Hij.component(a,b)*rij.getX(b);
            }
        }


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

    public void setAlpha(double alpha) { this.alpha = alpha; }

    public void setAvgs(double avgU, double avgF2, double avgH, double avg4) {
        this.avgU = avgU; this.avgF2 = avgF2; this.avgH = avgH; this.avg4 = avg4;}

    public void setUlat(double ULat) {
        this.Ulat = ULat;
    }
}

//double r6=r2*r2*r2;
//double w = -4/r6;
//double dw = 24/r6/r;
//
//double r02 = r0*r0;
//double r06=r02*r02*r02;
//
//double w0 = -4/r06;
//double dw0 = 24/r06/r0;
//
//double a = 0.1;
//double vr2 = -(Math.exp(beta*w) - a);
//double dvr2 = -beta*dw*Math.exp(beta*w);
//double vr20 = -(Math.exp(beta*w0) - a);
//double dvr20 = -beta*dw0*Math.exp(beta*w0);
//
//double vr2_fs = vr2 - vr20 + dvr20*(r0-r);
//double dvr2_fs = dvr2 - dvr20;
//
//sum += beta/2.0*vr2_fs/r2*Fij.dot(rij)/r + beta*dvr2_fs/r2;






// Jan 20



//        double L = box.getBoundary().getBoxSize().getX(0);
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            for (int j=0; j < dim; j++) {
//                double xi = atom.getPosition().getX(j);
//                double Fi = F0[iAtom].getX(j);
//                double Hii = array[dim*iAtom+j][dim*iAtom+j];
//
////                double n = 0.25;
////                double sinXi = Math.sin(n*2*Math.PI/L*xi);
////                double cosXi = Math.cos(n*2*Math.PI/L*xi);
////                double vi = sinXi;
////                double dvi = n*2*Math.PI/L*cosXi;
////                double vi = xi*xi;
////                double dvi = 2*xi;
////                java.util.Random r = new java.util.Random();
//
//                Random random = new Random();
//                double r = random.nextGaussian();
//                double vi = Fi + r;
//                double dvi = -Hii;
//                sumFv += Fi*vi;
//                sumH += dvi;
////                System.out.println(r);
//            }
//        }

//
//        Matrix www = eVecs.transpose().times(F_matrix);
//
//        double s =0;
//        for (int i=0; i<dof;i++) {
//            s += www.get(i, 0);
//        }
//        System.out.println("s: " + s);
//        for (int i = 0; i < numAtoms; i++) {
//            System.out.println(1/array[i*3][i*3] + "  " + F0[i].getX(0)/array[i*3][i*3]);
//        }

//        System.out.println("Hxx: " + array[0][0]);

//        int nn=0;
//        double sumFHiF = 0;
//        double[][] Hii_array = new double[dim][dim];
//        double[][] Fi_array = new double[dim][dim];
//        for (int k = 0; k < numAtoms; k++) {
//            for (int i = 0; i < dim; i++) {
//                Fi_array[i][0] = F0[k].getX(i);
//                for (int j = 0; j < dim; j++) {
//                    Hii_array[i][j] = array[k*dim+i][k*dim+j];
//                }
//            }
//            Matrix Hii = new Matrix(Hii_array);
//            Matrix Hii_inv = Hii.inverse();
//
//            Matrix Fi = new Matrix(Fi_array);
//            Matrix dot = Hii_inv.times(Fi);
//            Matrix dotdot = Fi.transpose().times(dot);
//            double tmp = Math.pow(90, 3);
//            if (Hii.det() > tmp) {
//                sumFHiF += dotdot.get(0,0);
//                nn++;
//            }
//        }

//        EigenvalueDecomposition ev_Hii_ = Hii.eig();
//        double[] ev_Hii = ev_Hii_.getRealEigenvalues();
////        System.out.println("ev: " + ev_Hii[0] + " " + ev_Hii[1] + " " + ev_Hii[2]);
//
//
////        System.out.println("dr: " + mapping.get(0,0) + " " + mapping.get(1,0) + " " + mapping.get(2,0));
//
//        System.out.println("Hxx: "+Hii.get(0,0) + "  " + Hii.get(1,1) + "  " + Hii.get(2,2));
//        System.out.println("dx: "+Fi.get(0,0)/Hii.get(0,0) + "  " + Fi.get(1,0)/Hii.get(1,1) + "  " + Fi.get(2,0)/Hii.get(2,2));
//        System.out.println();

//        System.out.println(array[n][n]   + " " + array[n][n+1]   + " " + array[n][n+2]);
//        System.out.println(array[n+1][n]   + " " + array[n+1][n+1]   + " " + array[n+1][n+2]);
//        System.out.println(array[n+2][n]   + " " + array[n+2][n+1]   + " " + array[n+2][n+2]);
//        System.out.println();
//double  test = 0;

// // 1D
//        int n = 0;
//        double sumFdr = 0;
//        double sumFdrHMA = 0;
//        double tmp = 0;
//        IAtomList atoms = box.getLeafList();
//        for (int i = 0; i < numAtoms; i++) {
//            IAtom atom = atoms.get(i);
//            for (int j = 0; j < dim; j++) {
//                double Fx0 = F0[i].getX(j);
//                double Hxx = array[dim*i+j][dim*i+j];
//                if (Hxx > tmp && step % 100 ==0) {
//                    double x0 = pos0[i].getX(j);
//                    double x1 = x0 + Fx0 / Hxx;
//                    Vector v = box.getSpace().makeVector();
//                    v.E(pos0[i]);
//                    v.setX(j, x1);
//                    atom.getPosition().E(v);
//                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
//                    atom.getPosition().PE(shift);
//                    pc.updateAtom(atom);
//                    pc.computeAll(true);
//                    double Fx1 = pc.getForces()[i].getX(j);
//                    double x2 = 0;
//                    double Fx2 = 0;
//
//                    double xTol = 0.00001;
//
//                    while (Math.abs(x1 - x0) > xTol) {
//                        double m = (Fx1 - Fx0) / (x1 - x0);
//                        x2 = x1 - Fx1 / m;
//                        v.E(pos0[i]);
//                        v.setX(j, x2);
//                        atom.getPosition().E(v);
//                        shift = box.getBoundary().centralImage(atom.getPosition());
//                        atom.getPosition().PE(shift);
//                        pc.updateAtom(atom);
//                        pc.computeAll(true);
//                        Fx2 = pc.getForces()[i].getX(j);
//                        x0 = x1;
//                        x1 = x2;
//                        Fx0 = Fx1;
//                        Fx1 = Fx2;
//                    }
//
//                    xmin[i].setX(j, x1);
//
//                    atom.getPosition().E(pos0[i]);
//                    pc.updateAtom(atom);
//                    pc.computeAll(false);
//                }
//
//                double dd = xmin[i].getX(j) - lattice[i].getX(j);
//
//                Vector dr = box.getSpace().makeVector();
//                dr.Ev1Mv2(pos0[i] , xmin[i]);
//                box.getBoundary().nearestImage(dr);
//                Vector drHMA = box.getSpace().makeVector();
//                drHMA.Ev1Mv2(pos0[i] , lattice[i]);
//                box.getBoundary().nearestImage(drHMA);
//
//                sumFdr += F0[i].getX(j)*dr.getX(j);
//                sumFdrHMA += F0[i].getX(j)*drHMA.getX(j);
//                n++;
//            }
//        }

//test = 0;
//FD
//        for (int i = 0; i < numAtoms; i++) {
//            for (int j = 0; j < dim; j++) {
//                double Hxx = array[dim*i+j][dim*i+j];
//                if (Hxx > tmp) {
//                    sumF2 += F0[i].getX(j)*F0[i].getX(j);
//                    sumHxx += Hxx;
//                    sF2_H += F0[i].getX(j)*F0[i].getX(j)/Hxx;
//                    n++;
//                }
//            }
//        }
//        int i = 0;
//        int j = 0;
//        double ddx = lattice[i].getX(j) - box.getLeafList().get(i).getPosition().getX(j) - F0[i].getX(j)/array[dim*i+j][dim*i+j];
//        System.out.println(ddx);
//
//        double dx_max = 0.5;
//        int len = 50;
//        int iAtom = 0;
//        IAtom atom = box.getLeafList().get(0);
//        double Hxx = array[dim*iAtom][dim*iAtom];
//
//        if (Hxx > 0 && step%10 == 0) {
//            for (int i = -len; i <= len; i++) {
//                Vector dx = box.getSpace().makeVector();
//                dx.setX(0, i*dx_max/len);
//                atom.getPosition().Ev1Pv2(pos0[iAtom], dx);
//                pc.updateAtom(atom);
//                double Ui = pc.computeAll(false);
//                System.out.println(dx.getX(0) + " " + (Ui-U0)/numAtoms + "      " + Hxx);
//            }
//            atom.getPosition().E(pos0[iAtom]);
//            pc.updateAtom(atom);
//            pc.computeAll(false);
//
//            System.out.println("");
//
//        }
//
//        double d = 1e-3;
//        Vector Fpx = box.getSpace().makeVector();
//        Vector Fmx = box.getSpace().makeVector();
//        double sumHxxx = 0;
//        for (IAtom atom : box.getLeafList()) {
//            for (int j = 0; j < dim; j++) {
//                Vector dx = box.getSpace().makeVector();
//                dx.setX(j, d);
//                int iAtom = atom.getLeafIndex();
//                atom.getPosition().PE(dx);
//                pc.updateAtom(atom);
//                pc.computeAll(true);
//                Fpx.E(pc.getForces()[iAtom]);
//
//                atom.getPosition().PEa1Tv1(-2, dx);
//                pc.updateAtom(atom);
//                pc.computeAll(true);
//                Fmx.E(pc.getForces()[iAtom]);
//
//                double Hxx  = array[dim*iAtom+j][dim*iAtom+j];
//                double Hxxx = -(Fpx.getX(j) -2*F0[iAtom].getX(j) + Fmx.getX(j))/(d*d);
//                if (Hxx > tmp) {
//                    sumHxxx += F0[iAtom].getX(j) * Hxxx / Hxx / Hxx;
//                }
//                atom.getPosition().E(pos0[iAtom]);
//                pc.updateAtom(atom);
//                pc.computeAll(false);
//            }
//        }


//3D

//3D
//test = 0;

//        Vector sumdr = box.getSpace().makeVector();
//        int n = 0;
//        double sumFdr = 0;
//        double sumFdrHMA = 0;
//        double tmp = 0;
//        double fac = 0.001;
//
//        IAtomList atoms = box.getLeafList();
//        Vector dr = box.getSpace().makeVector();
//        Vector rOld = box.getSpace().makeVector();
//        Vector fOld = box.getSpace().makeVector();
//        for (int i = 0; i < numAtoms; i++) {
//            double[][] Hii_array = new double[dim][dim];
//            for (int a = 0; a < dim; a++) {
//                for (int b = 0; b < dim; b++) {
//                    Hii_array[a][b] = array[i * dim + a][i * dim + b];
//                }
//            }
//            Matrix Hii = new Matrix(Hii_array);
//            EigenvalueDecomposition ed = Hii.eig();
//            double[] eVals = ed.getRealEigenvalues();
//
//            if (eVals[0] < tmp || eVals[1] < tmp || eVals[2] < tmp) continue;
//
//            if (step % 20 == 0) {
//                IAtom atom = atoms.get(i);
//                rOld.E(pos0[i]);
//                fOld.E(F0[i]);
//                for (int k = 0; k < 20; k++) {
//                    dr.Ea1Tv1(fac, fOld);
//                    atom.getPosition().Ev1Pv2(rOld, dr);
//                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
//                    atom.getPosition().PE(shift);
//                    pc.updateAtom(atom);
//                    pc.computeAll(true);
//                    fOld.E(pc.getForces()[i]);
//                    rOld.E(atom.getPosition());
//                }
//
//                rmin[i].E(rOld);
//
//                atom.getPosition().E(pos0[i]);
//                pc.updateAtom(atom);
//                pc.computeAll(false);
//                n++;
//            }
//
//            Vector rrmin = box.getSpace().makeVector();
//            Vector rrminHMA = box.getSpace().makeVector();
//            rrmin.Ev1Mv2(pos0[i], rmin[i]);
//            rrminHMA.Ev1Mv2(pos0[i], lattice[i]);
//            box.getBoundary().nearestImage(rrmin);
//            box.getBoundary().nearestImage(rrminHMA);
//
//            sumFdr += F0[i].dot(rrmin);
//            sumFdrHMA += F0[i].dot(rrminHMA);
//            n++;
//        }







//dN
//        String method = "RMSprop";
//        String method = "SD";
//        String method = "XXX";
//
//        int maxIter = 1000;
//        double totalD = 0;
//        double Umin = 0;
//        if (method == "RMSprop") {
//            double eta = 0.01;
//            double eps = 1e-8;
//
//            double[] Gt = new double[dim * (numAtoms - 1)];
//            for (int iter = 0; iter < maxIter; iter++) {
//                totalD = 0;
//                double sumDir2 = 0;
//                for (int i = 1; i < numAtoms; i++) {
//                    Vector dr = box.getSpace().makeVector();
//                    for (int j = 0; j < dim; j++) {
//                        int ij = dim * (i - 1) + j;
//                        double gt = pc.getForces()[i].getX(j);
//                        Gt[ij] += gt * gt;
////                        double dir = eta * gt;
//                        double dir = eta / (eps + Math.sqrt(Gt[ij])) * gt;
//                        dr.setX(j, dir);
//                        sumDir2 += dir*dir;
//                        totalD -= gt*dir;
//                    }//j
//                    atoms.get(i).getPosition().PE(dr);
//                    Vector shift = box.getBoundary().centralImage(atoms.get(i).getPosition());
//                    atoms.get(i).getPosition().PE(shift);
//                }//i
//                totalD /= Math.sqrt(sumDir2);
////                System.out.println((iter+1) + "  " + pc.getLastEnergy() + "  " + totalD);
//
//
//                pc.init();
//                pc.computeAll(true);
//            }//t
//
//
//        } else if (method == "SD") {
//
//
//            SteepestDescent sd = new SteepestDescent(new FunctionMultiDimensionalDifferentiable() {
//                @Override
//                public double df(int[] d, double[] x) {
//                    throw new RuntimeException("This is not the method you're looking for");
//                }
//
//                protected void assignCoords(double[] x) {
//                    for (int i = 1; i < numAtoms; i++) {
//                        for (int j = 0; j < dim; j++) {
//                            atoms.get(i).getPosition().setX(j, x[dim * (i - 1) + j]);
//                        }
//                        Vector shift = box.getBoundary().centralImage(atoms.get(i).getPosition());
//                        atoms.get(i).getPosition().PE(shift);
//                    }
//                    pc.init();
//                }
//
//                @Override
//                public double[] gradf(double[] x) {
//                    assignCoords(x);
//                    pc.computeAll(true);
//                    Vector[] forces = pc.getForces();
//                    double[] rv = new double[x.length];
//                    for (int i = 1; i < numAtoms; i++) {
//                        for (int j = 0; j < dim; j++) {
//                            rv[dim * (i - 1) + j] = -forces[i].getX(j);
//                        }
//                    }
//                    return rv;
//                }
//
//                @Override
//                public double f(double[] x) {
//                    assignCoords(x);
//                    double u = pc.computeAll(false);
//                    return u;
//                }
//
//                @Override
//                public int getDimension() {
//                    return dim * (numAtoms - 1);
//                }
//            });
//
//            double[] x0 = new double[dim * (numAtoms - 1)];
//            double[] xStep = new double[dim * (numAtoms - 1)];
//            for (int i = 1; i < numAtoms; i++) {
//                for (int j = 0; j < dim; j++) {
//                    x0[dim * (i - 1) + j] = atoms.get(i).getPosition().getX(j);
//                    xStep[dim * (i - 1) + j] = 0.01;
//                }
//            }
//            if (pc != null) {
//                double tol = 1e-4;
//                sd.minimize(x0, xStep, tol, maxIter);
//                totalD = sd.getTotalD();
//            }
//
//            Umin = pc.getLastEnergy();
//            double T = 1/beta;
//            System.out.println(step + " " + Umin/numAtoms + "    " + totalD);
//
//            for (IAtom atom : atoms) {
//                int iAtom = atom.getLeafIndex();
//                rmin[iAtom].E(atom.getPosition());
//                atom.getPosition().E(pos0[iAtom]);
//            }
//            pc.init();
//            pc.computeAll(false);
//
//        }



// COM




//        double sumFdrHMA = 0;
//        double sumFx2 = 0;
//        double sumxHx = 0;
//
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            Vector drHMA = box.getSpace().makeVector();
//            drHMA.Ev1Mv2(pos0[iAtom], rmin[iAtom]);
//            drHMA.ME(com);
//            box.getBoundary().nearestImage(drHMA);
//
//            sumFdrHMA += F0[iAtom].dot(drHMA);
//            for (int a = 0; a < dim; a++) {
//                int ij = iAtom*dim+a;
//                sumxHx += drHMA.getX(a)*array[ij][ij]*drHMA.getX(a);
//                sumFx2 += Math.pow(F0[iAtom].getX(a)*drHMA.getX(a), 2);
//            }
//
//        }
//
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            Vector drHMAhalf = box.getSpace().makeVector();
//            drHMAhalf.Ev1Mv2(pos0[iAtom], rmin[iAtom]);
//            drHMAhalf.ME(com);
//            box.getBoundary().nearestImage(drHMAhalf);
//            drHMAhalf.TE(0.5);
//            atom.getPosition().Ev1Pv2(rmin[iAtom], drHMAhalf);
//        }
//        pc.init();
//        pc.computeAll(true);
//        Vector[] Fhalf = pc.getForces();
//
//        double sumFdrHMAhalf = 0;
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            Vector drHMAhalf = box.getSpace().makeVector();
//            drHMAhalf.Ev1Mv2(atom.getPosition(), rmin[iAtom]);
//            drHMAhalf.ME(com);
//            box.getBoundary().nearestImage(drHMAhalf);
//            sumFdrHMAhalf += Fhalf[iAtom].dot(drHMAhalf);
//        }
//
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            rmin[iAtom].E(atom.getPosition());
//            atom.getPosition().E(pos0[iAtom]);
//        }
//        pc.init();
//        pc.computeAll(false);
//
//
//        System.out.println(step + "  " + beta*sumFdrHMA/(dim*(numAtoms-1)) + "  " + beta*sumFdrHMAhalf/(dim*(numAtoms-1)));


//        double y1 = sumFdrHMAhalf/numAtoms;
//
//        double alpha_ = 0.15;
//        double y2 = y1 - alpha_*(sumFdrHMA + dim*(numAtoms-1.0)/beta)/numAtoms;

//        System.out.println(step + " " + y1 + " " + y2 + "     " + sumFdrHMA/numAtoms + "  " + totalD);

//        x[0] = U0/numAtoms;

//        double sumdUi = 0;
//        double dU = 0;
//        if (step > 0) {
//            for (IAtom atom : atoms) {
//                int iAtom = atom.getLeafIndex();
//                Vector dx = box.getSpace().makeVector();
//                dx.Ev1Mv2(atom.getPosition(), rOld[iAtom]);
//                box.getBoundary().nearestImage(dx);
//                sumdUi -= F0[iAtom].dot(dx) + Fold[iAtom].dot(dx);
//            }
//            sumdUi /= (2.0*numAtoms);
//            dU = (U0-Uold)/numAtoms;
//            double y0 = Uold/numAtoms;
//            double fac = -0.31338696366356233;
//            System.out.println(step + "  " + y0 + "  " + sumdUi + "  " + dU);
////            System.out.println(step + "  " + y0 + "  " + (y0-fac*sumdUi));
//
//        }

//        System.out.println(step + "  " + dU + "  " + sumdUi + "  " + (dU-sumdUi));

//        x[1] = dim*(numAtoms-1.0)/numAtoms/2.0/beta + U0/numAtoms + 1.0/2.0*sumFdrHMA/numAtoms;
//        x[2] = dim*(numAtoms-1.0)/numAtoms/6.0/beta - 4.0/3.0*y2 + U0/numAtoms + 1.0/6.0*sumFdrHMA/numAtoms + 4.0/3.0*sumFdrHMAhalf/numAtoms;
//        x[3] = sumFdrHMAhalf/numAtoms;
//        x[4] = sumFdrHMA/numAtoms;

//        double y3 = dim*(numAtoms-1.0)/numAtoms/6.0/beta + U0/numAtoms + 1.0/6.0*sumFdrHMA/numAtoms + 4.0/3.0*sumFdrHMAhalf/numAtoms;
//        System.out.println(step + "  " + x[0] + "   " + x[2] + "   " + y3 + "   " + totalD);


//        double sumFdrMin = 0, sumFdr = 0;
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//
//            Vector dr = box.getSpace().makeVector();
//            dr.Ev1Mv2(atom.getPosition(), lattice[iAtom]);
//            dr.ME(com);
//            box.getBoundary().nearestImage(dr);
//
//            Vector drMin = box.getSpace().makeVector();
//            drMin.Ev1Mv2(atom.getPosition(), rmin[iAtom]);
//            drMin.ME(com);
//            box.getBoundary().nearestImage(drMin);
//
//            sumFdr += F0[iAtom].dot(dr);
////            sumFdrMin += F0[iAtom].dot(drMin);
//            sumFdrMin += F0[iAtom].dot(atom.getPosition());
//
//        }

//        System.out.println(step + "  " + U0/numAtoms + "  " + sumFdrMin/numAtoms);

//        double uConv = U0/numAtoms;
//        double uHMA = dim*(numAtoms-1.0)/2.0/beta/numAtoms + U0/numAtoms + 1.0/2.0*sumFdr/numAtoms;
//        System.out.println(step + "  " + uConv + "  " + uHMA);

//        maxIter = 1000;
//        double hn=0.001;
//        double sumg2 = 0;
//        if (step % 1 == 0) {
//            double Uold = pc.getLastEnergy();
//            for (int t = 1; t <= maxIter; t++) {
//                Vector[] post = box.getSpace().makeVectorArray(numAtoms);
//                for (int i = 0; i < numAtoms; i++) {
//                    post[i].E(atoms.get(i).getPosition());
//                }
//                double maxFi2=-1.0;
//                for (int i = 0; i < numAtoms; i++) {
//                    if (pc.getForces()[i].squared() > maxFi2) {
//                        maxFi2 = pc.getForces()[i].squared();
//                    }
//                }
//
//                for (int i = 0; i < numAtoms; i++) {
//                    Vector dr = box.getSpace().makeVector();
//                    dr.Ea1Tv1(hn/Math.sqrt(maxFi2), pc.getForces()[i]);
//                    atoms.get(i).getPosition().PE(dr);
//                    Vector shift = box.getBoundary().centralImage(atoms.get(i).getPosition());
//                    atoms.get(i).getPosition().PE(shift);
//                }//i
//
//                pc.init();
//                double Unew = pc.computeAll(false);
//                if (Unew < Uold) {
//                    hn *= 1.2;
//                } else {
//                    hn *= 0.2;
//                    for (IAtom atom : atoms) {
//                        atom.getPosition().E(post[atom.getLeafIndex()]);
//                    }
//                    pc.init();
//                }
//                pc.computeAll(true);
//            }//t
//
//            for (IAtom atom : atoms) {
//                int iAtom = atom.getLeafIndex();
//                sumg2 += pc.getForces()[iAtom].squared();
//                rmin[iAtom].E(atom.getPosition());
//                atom.getPosition().E(pos0[iAtom]);
//            }
//            pc.init();
//            pc.computeAll(false);
//        } //if
//
//
//        double sumFdr = 0;
//        double sumFdrHMA = 0;
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            Vector dr = box.getSpace().makeVector();
//            Vector drHMA = box.getSpace().makeVector();
//            dr.Ev1Mv2(pos0[iAtom], rmin[iAtom]);
//            drHMA.Ev1Mv2(pos0[iAtom], lattice[iAtom]);
//            box.getBoundary().nearestImage(dr);
//            box.getBoundary().nearestImage(drHMA);
//            sumFdr += F0[iAtom].dot(dr);
//            sumFdrHMA += F0[iAtom].dot(drHMA);
//        }

//        double scale = n*1.0/(numAtoms);
//



//        x[1] = 1.0/2.0*sumFHiF/numAtoms - avgF2;
//        x[2] = U0/numAtoms + scale_nn*dim/2.0/beta - 1.0/2.0*sumFHiF/numAtoms;

//        if (avgU != 0) {
//            System.out.println(step + "  " + U0/numAtoms + " " + x[2]);
//        }

//        System.out.println(sumFx2/sumHxx + "  " + sumFy2/sumHyy + "  " + sumFz2/sumHzz);
//        System.out.println((1.0/beta - sF2_H));
//        System.out.println(step + "  "+ U0/numAtoms + "  " + sF2_H/(dof));
//        x[2] = sF2_H/(dof);
//        System.out.println(3.0/2.0*x[2]);

//        x[0] = U0/numAtoms; // <U/N>
//        x[1] = sumF2/numAtoms; // <F2/N>
//        x[2] = sumHxx/numAtoms; // <trH/N>
//        x[3] = sumHxx/sumF2; //<trH/F2>


//        this.avgU = avgU; this.avgF2 = avgF2; this.avgH = avgH; this.avg4 = avg4;}

//        x[0] = U0/numAtoms - avgU;
//        x[1] = sumF2/numAtoms - avgF2;
//        x[2] = sumHxx/numAtoms - avgH;
//        x[3] = sumHxx/sumF2 - avg4;

//        System.out.println(step + "  " + x[0] + "  " + x[1]);


//        x[0] = U0/numAtoms;
//        x[1] = (U0 - 1.0/2.0*F0[0].getX(0)*F0[0].getX(0)/array[0][0])/numAtoms;
//        x[1] = (U0 - 1.0/2.0/beta*sumHxx/Hxx_avg)/numAtoms;
//        x[1] = (U0 - 1.0/2.0/beta*sumHxy/Hxx_avg)/numAtoms;
//        x[1] = (U0 - 1.0/2.0*dotdot.get(0,0))/numAtoms;


//        System.out.println(step + " " + x[0]+ " " + x[1]+ " " + x[2]);
//        System.out.println(step + " " + U0 + "  " + beta*sumF2/(dof) + "  " + sumHxx/(dof));
//        System.out.println(step + " " + x[0] + "  " + x[1]);
//        step++;
//        double sumLambda = 0;
//        double sumLambda2 = 0;
//        for (int i=0; i<dof; i++) {
//            sumLambda += eVals[i];
//            sumLambda2 += eVals[i]*eVals[i];
//        }

//        double sumFHF = 0;
//        for (int i = 0; i < numAtoms; i++){
//            for (int j = 0; j < numAtoms; j++){
//                for (int a = 0; a < dim; a++) {
//                    for (int b = 0; b < dim; b++) {
//                        sumFHF += F0[i].getX(a)*array[dim*i+a][dim*j+b]*F0[j].getX(b);
//                    }
//                }
//            }
//        }

//        x[5] = (U0 + alpha*(sumLambda-beta*sumF2))/numAtoms;
//        double Hxx_avg = 322.0;
//        alpha = 1.0/2.0/beta/Hxx_avg;
//        x[1] = (dof/2.0/beta + U0 - alpha*beta*sumF2)/numAtoms;


//        x[1] = (U0 + alpha*(sumLambda-beta*sumF2))/numAtoms;
//        x[2] = sumHxx/dof;

//        x[0] = sumF2;
//        x[1] = sumLambda;
//        x[2] = sumLambda2;
//        x[3] = sumFHF;

//        x[4] = U0/numAtoms;
////        x[5] = (U0 + alpha*(sumLambda-beta*sumF2))/numAtoms;
//        double Hxx_avg = 3.357016516983524e+02;
//        alpha = 1.0/2.0/beta/Hxx_avg;
//        x[5] = (dof/2.0/beta + U0 - alpha*sumF2)/numAtoms;

//
//        for (IAtom atom : box.getLeafList()) {
//            atom.getPosition().PEa1Tv1(alpha, F0[atom.getLeafIndex()]);
//            Vector shift = box.getBoundary().centralImage(atom.getPosition());
//            atom.getPosition().PE(shift);
//        }
//        pc.init();
//        pc.computeAll(true);
//        double U = pc.getLastEnergy();
//        Vector[] F = pc.getForces();
//
//        double bU = beta1*U;
//        double sumFF0 = 0;
//        for (IAtom atom : box.getLeafList()) {
//            sumFF0 += F0[atom.getLeafIndex()].dot(F[atom.getLeafIndex()]);
//        }
//
//        double sumL = 0;
//        for (int i = 0; i < dof; i++) {
//            sumL += eVals[i]/(1-alpha*eVals[i]);
//        }
//
//        x[4] = (beta1-beta)*U0/numAtoms; // du
//        x[5] = (bU - b0U - Math.log(J))/numAtoms; // phi
//        x[6] = (-beta1*sumFF0 + sumL)/numAtoms; //dphi/dalpha
//        x[7] = U0/numAtoms; //u0
//
//        for (IAtom atom : box.getLeafList()) {
//            atom.getPosition().E(pos0[atom.getLeafIndex()]);
//        }
//
//        pc.init();
//        pc.computeAll(true);
//
//        Vector dr = box.getSpace().makeVector();
//        double scale = Math.sqrt(beta/beta1);
//        for (IAtom atom : box.getLeafList()) {
//            dr.Ev1Mv2(atom.getPosition(), lattice[atom.getLeafIndex()]);
//            box.getBoundary().nearestImage(dr);
//            atom.getPosition().E(lattice[atom.getLeafIndex()]);
//            atom.getPosition().PEa1Tv1(scale, dr);
//            Vector shift = box.getBoundary().centralImage(atom.getPosition());
//            atom.getPosition().PE(shift);
//        }
//
//        pc.init();
//        pc.computeAll(false);
//        double U2 = pc.getLastEnergy();
//        double bU2 = beta1*U2;
//        J = Math.pow(beta/beta1, dim*(numAtoms-1.0)/2.0);
//
//        x[8] = (bU2 - b0U - Math.log(J))/numAtoms; // phi
//
//        for (IAtom atom : box.getLeafList()) {
//            atom.getPosition().E(pos0[atom.getLeafIndex()]);
//        }
//        pc.init();
//        pc.computeAll(true);
//
//        for (IAtom atom : atoms) {
//            int iAtom = atom.getLeafIndex();
//            rOld[iAtom].E(atom.getPosition());
//            Fold[iAtom].E(pc.getForces()[iAtom]);
//        }
//        Uold = pc.getLastEnergy();


//        EigenvalueDecomposition ed = matrix.eig();
//        double[] eVals = ed.getRealEigenvalues();




/*** bU-b0U0=0 **/
//double[][] arrayFF = new double[dof][dof];
//double[][] arrayFFH = new double[dof][dof];
//double[][] arrayJ = new double[dof][dof];

//double prodJ = 1;
//double sumHii = 0;
//double sumHij = 0;
//
//alpha = 0.00001*dbeta;
//
//        for (int i=0; i<dof; i++) {
//prodJ *= (1-alpha*eVals[i]);
//sumHii += arrayH[i][i];
//        for (int j=0; j<dof; j++) {
//sumHij += arrayH[i][j]*arrayH[i][j];
//        }
//        }
//
//
//double a = -alpha*sumHii;
//double b = -alpha*sumHii-0.5*alpha*alpha*sumHij;
//double U0 = pc.getLastEnergy();
//double bU0 = beta*U0;
//Vector[] F0 = box.getSpace().makeVectorArray(numAtoms);
//double sumF2 = 0;
//        for (int i=0; i<numAtoms; i++) {
//F0[i].E(pc.getForces()[i]);
//sumF2 += F0[i].squared();
//        }


//EigenvalueDecomposition ed = matrix.eig();
//double[] eVals = ed.getRealEigenvalues();

package etomica.normalmode;

import etomica.atom.Atom;
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
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import org.apache.xml.utils.IntVector;


public class MeterFM implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pc;
    protected double beta;
    protected double trH;
    protected double[][] array;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected int numAtoms;
    protected int step = 0;
    protected double alpha;
    protected Vector[] pos0, lattice, xmin, rmin;
    protected double dbeta;
    protected double avgU, avgF2, avgH, avg4;

    public MeterFM(PotentialCompute pc, double temperature, Box box, double alpha, double dbeta) {
        int nData = 3;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("FM",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pc = pc;
        this.beta = 1/temperature;
        this.box = box;
        this.alpha = alpha;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        pos0 = box.getSpace().makeVectorArray(numAtoms);
        lattice = box.getSpace().makeVectorArray(numAtoms);
        xmin = box.getSpace().makeVectorArray(numAtoms);
        rmin = box.getSpace().makeVectorArray(numAtoms);

        for (IAtom atom : box.getLeafList()) {
            lattice[atom.getLeafIndex()].E(atom.getPosition());
        }
        this.dbeta = dbeta;

        this.avgU = 0;
        this.avgF2 = 0;
        this.avgH = 0;
        this.avg4 = 0;
    }

    @Override
    public IData getData() {
        double beta1 = beta + dbeta;
        double[] x = data.getData();

        for (IAtom atom : box.getLeafList()) {
            pos0[atom.getLeafIndex()].E(atom.getPosition());
        }
        trH = 0;
        array = new double[dim*numAtoms][dim*numAtoms];
        pc.computeAll(true, this);
        double U0 = pc.getLastEnergy();
        double b0U = beta*U0;
        Vector[] F0 = box.getSpace().makeVectorArray(numAtoms);
        for (int i=0; i<numAtoms; i++) {
            F0[i].E(pc.getForces()[i]);
        }

//        Matrix matrix = new Matrix(array);
//        EigenvalueDecomposition ed = matrix.eig();
//        double[] eVals = ed.getRealEigenvalues();
//        Matrix eVecs = ed.getV();
//
//        double J = 1;
//        for (int i = 0; i < array.length; i++) {
//            J *= (1.0 - alpha*eVals[i]);
//        }
//
//        Matrix www = eVecs.transpose().times(F_matrix);
//
//        double s =0;
//        for (int i=0; i<array.length;i++) {
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
        double  test = 0;

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

        test = 0;
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
        test = 0;

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
        int iter = 1000;
        double eps = 1e-8;
        double eta = 0.01;
        IAtomList atoms = box.getLeafList();

        double sumg2 = 0;
        if (step % 1 == 0) {
            double[] Gt = new double[dim*numAtoms];
            for (int t = 1; t <= iter; t++) {
                pc.computeAll(true);
                for (int i = 0; i < numAtoms; i++) {
                    Vector dr = box.getSpace().makeVector();
                    for (int j = 0; j < dim; j++) {
                        int ij = dim * i + j;
                        double gt = pc.getForces()[i].getX(j);
                        Gt[ij] += gt * gt;
                        dr.setX(j, eta / Math.sqrt(eps + Gt[ij]) * gt);
                        if (t == iter) sumg2 += gt * gt;
                    }
                    atoms.get(i).getPosition().PE(dr);
                    Vector shift = box.getBoundary().centralImage(atoms.get(i).getPosition());
                    atoms.get(i).getPosition().PE(shift);
                }
                pc.init();
            } //t

            for (IAtom atom : atoms) {
                int iAtom = atom.getLeafIndex();
                rmin[iAtom].E(atom.getPosition());
                atom.getPosition().E(pos0[iAtom]);
            }
            pc.init();
            pc.computeAll(false);
        } //if

        double sumFdr = 0;
        double sumFdrHMA = 0;
        for (IAtom atom : atoms) {
            int iAtom = atom.getLeafIndex();
            Vector dr = box.getSpace().makeVector();
            Vector drHMA = box.getSpace().makeVector();
            dr.Ev1Mv2(pos0[iAtom], rmin[iAtom]);
            drHMA.Ev1Mv2(pos0[iAtom], lattice[iAtom]);
            box.getBoundary().nearestImage(dr);
            box.getBoundary().nearestImage(drHMA);
            sumFdr += F0[iAtom].dot(dr);
            sumFdrHMA += F0[iAtom].dot(drHMA);
        }

//        double scale = n*1.0/(numAtoms);
//
        x[0] = U0/numAtoms;
        x[1] = dim/2.0/beta + U0/numAtoms + 1.0/2.0*sumFdr/numAtoms;
        x[2] = dim/2.0/beta + U0/numAtoms + 1.0/2.0*sumFdrHMA/numAtoms;


//        System.out.println(step + "   " + x[0] + "  " + x[1] + "    sum_g2: " + sumg2);
//        System.out.println(step + "   " + x[0] + "  " + x[1] + " " + x[2]);

        System.out.println(step + "   " + x[0] + "  " + x[1] + " " + x[2] + " sum_g2: " + sumg2);




//        x[1] = 1.0/2.0*sumFHiF/numAtoms - avgF2;
//        x[2] = U0/numAtoms + scale_nn*dim/2.0/beta - 1.0/2.0*sumFHiF/numAtoms;

//        if (avgU != 0) {
//            System.out.println(step + "  " + U0/numAtoms + " " + x[2]);
//        }

//        System.out.println(sumFx2/sumHxx + "  " + sumFy2/sumHyy + "  " + sumFz2/sumHzz);
//        System.out.println((1.0/beta - sF2_H));
//        System.out.println(step + "  "+ U0/numAtoms + "  " + sF2_H/(dim*numAtoms));
//        x[2] = sF2_H/(dim*numAtoms);
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
//        System.out.println(step + " " + U0 + "  " + beta*sumF2/(dim*numAtoms) + "  " + sumHxx/(dim*numAtoms));
//        System.out.println(step + " " + x[0] + "  " + x[1]);
//        step++;
//        double sumLambda = 0;
//        double sumLambda2 = 0;
//        for (int i=0; i<array.length; i++) {
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
//        x[1] = (dim*numAtoms/2.0/beta + U0 - alpha*beta*sumF2)/numAtoms;


//        x[1] = (U0 + alpha*(sumLambda-beta*sumF2))/numAtoms;
//        x[2] = sumHxx/array.length;

//        x[0] = sumF2;
//        x[1] = sumLambda;
//        x[2] = sumLambda2;
//        x[3] = sumFHF;

//        x[4] = U0/numAtoms;
////        x[5] = (U0 + alpha*(sumLambda-beta*sumF2))/numAtoms;
//        double Hxx_avg = 3.357016516983524e+02;
//        alpha = 1.0/2.0/beta/Hxx_avg;
//        x[5] = (dim*numAtoms/2.0/beta + U0 - alpha*sumF2)/numAtoms;

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
//        for (int i = 0; i < array.length; i++) {
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

        step++;
        return data;
    }

    public void fieldComputeHessian(int i, Tensor Hii) { }

    public void pairComputeHessian(int i, int j, Tensor Hij) {
        trH += 2*Hij.trace();
        for (int a = 0; a < dim; a++) {
            for (int b = 0; b < dim; b++) {
                array[dim*i+a][dim*j+b] = -Hij.component(a,b);
                array[dim*j+b][dim*i+a] = -Hij.component(a,b);
                //self
                array[dim*i+a][dim*i+b] += Hij.component(a,b);
                array[dim*j+a][dim*j+b] += Hij.component(a,b);
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
}
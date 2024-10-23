package etomica.normalmode;

import etomica.atom.IAtom;
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
    protected double alpha;
    protected Vector[] pos0;

    public MeterFM(PotentialCompute pc, double temperature, Box box, double alpha) {
        int nData = 8;
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

    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        for (IAtom atom : box.getLeafList()) {
            pos0[atom.getLeafIndex()].E(atom.getPosition());
        }
        trH = 0;
        array = new double[dim*numAtoms][dim*numAtoms];
        pc.computeAll(true, this);
        double U0 = pc.getLastEnergy();
        double b0U = beta*U0;
        Vector[] F0 = pc.getForces();
        Matrix matrix = new Matrix(array);
        EigenvalueDecomposition ed = matrix.eig();
        double[] eVals = ed.getRealEigenvalues();
        double J = 1;
        for (int i = 0; i < array.length; i++) {
            J *= (1.0 - alpha*eVals[i]);
        }

        double sumF2 = 0;
        for (int i=0; i<numAtoms; i++) {
            sumF2 += F0[i].dot(F0[i]);
        }
        double sumLambda = 0;
        double sumLambda2 = 0;
        for (int i=0; i<array.length; i++) {
            sumLambda += eVals[i];
            sumLambda2 += eVals[i]*eVals[i];
        }
        double sumFHF = 0;
        for (int i = 0; i < numAtoms; i++){
            for (int j = 0; j < numAtoms; j++){
                for (int a = 0; a < dim; a++) {
                    for (int b = 0; b < dim; b++) {
                        sumFHF += F0[i].getX(a)*array[dim*i+a][dim*j+b]*F0[j].getX(b);
                    }
                }
            }
        }
        x[0] = sumF2;
        x[1] = sumLambda;
        x[2] = sumLambda2;
        x[3] = sumFHF;

        for (IAtom atom : box.getLeafList()) {
            atom.getPosition().PEa1Tv1(alpha, F0[atom.getLeafIndex()]);
        }
        pc.init();

        pc.computeAll(true);
        double U = pc.getLastEnergy();
        Vector[] F = pc.getForces();
        double beta1 = beta + 0.01;
        double bU = beta1*U;
        double sumFF0 = 0;
        for (IAtom atom : box.getLeafList()) {
            sumFF0 += F0[atom.getLeafIndex()].dot(F[atom.getLeafIndex()]);
        }

        double sumL = 0;
        for (int i = 0; i < array.length; i++) {
            sumL += eVals[i]/(1-alpha*eVals[i]);
        }

        x[4] = (beta1-beta)*U0/numAtoms;
        x[5] = (bU - b0U - Math.log(J))/numAtoms;
        x[6] = (-beta1*sumFF0 + sumL)/numAtoms;
        x[7] = U0/numAtoms;

//        System.out.println(step + " " + (beta1-beta)*U0+ "   "+ x[0] + "   J: " + Math.pow(J, 1.0/(dim*numAtoms)));
//        step++;
//        System.out.println((bU-b0U) + " " + -Math.log(J) + "  " + x[0]);
//        System.out.println(x[0] + "  J: " + J);
//        System.out.println(x[0] + "  J: " + J);

        for (IAtom atom : box.getLeafList()) {
            atom.getPosition().E(pos0[atom.getLeafIndex()]);
        }
        pc.init();

//        double sumF2 = 0;
//        for (IAtom atom : box.getLeafList()) {
//            sumF2 += F0[atom.getLeafIndex()].dot(F0[atom.getLeafIndex()]);
//        }
//        x[0] = U0/numAtoms;
//        x[1] = (U0 + alpha*(trH - beta*sumF2))/numAtoms;
//        x[2] = (trH - beta*sumF2)/numAtoms;
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
}

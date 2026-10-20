/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;

/**
 * Computes various quantities that ZENO computes
 * "eq 21" refers to eq 21 from <a href="https://dx.doi.org/10.1103/PhysRevE.78.046712">Improved path integration
 * method for estimating the intrinsic viscosity of arbitrarily shaped particles</a>
 */
public class MeterZENO implements IAction, IDataSource {
    protected final Box box;
    public WalkerExterior walker;
    protected long numWalks;
    protected IRandom random;
    protected Vector KPlus;
    protected Vector KMinus;
    protected Tensor VPlus;
    protected Tensor VMinus;
    protected Vector uAvg;
    protected Tensor qvAvg;
    protected Tensor qwAvg;
    protected double tAvg;
    protected long hits;
    protected long totalWalks, totalSteps;
    double[] sigmaByType;
    protected boolean eigenvaluesForPade;

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final boolean doEq21;
    protected WalkerExterior.Shape shape;

    public MeterZENO(Box box, IRandom random, double[] sigmaByType, boolean doEq21, WalkerExterior.Shape shape) {
        this.box = box;
        this.random = random;
        BoundingSphereGenerator.BoundingSphere bs  = BoundingSphereGenerator.getBoundingSphere(box, sigmaByType, shape);
        Vector boundingSphereCenter = bs.center;
        double boundingSphereRadius = bs.radius;

        this.walker = new WalkerExterior(box, sigmaByType, random, boundingSphereRadius, boundingSphereCenter, shape);
        this.KPlus = new Vector3D();
        this.KMinus = new Vector3D();
        this.VPlus = new Tensor3D();
        this.VMinus = new Tensor3D();
        this.sigmaByType = sigmaByType;

        // hydrodynamic radius, intrinsic conductivity, practical conductivity, q_eta, intrinsic viscosity, practical viscosity, practical viscosity just W
        data = new DataDoubleArray(7);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("zeno", Null.DIMENSION, new int[]{7});
        tag = new DataTag();
        dataInfo.addTag(tag);

        this.doEq21 = doEq21;
        this.shape = shape;
        if (doEq21) {
            uAvg = new Vector3D();
            qvAvg = new Tensor3D();
            qwAvg = new Tensor3D();
        }
    }

    public void setEigenvaluesForPade(boolean b) {
        eigenvaluesForPade = b;
    }

    public void setBoundingSphere(Vector center, double radius) {
        walker.setBoundingSphere(center, radius);
    }

    public void resetBoundingSphere() {
        BoundingSphereGenerator.BoundingSphere bs  = BoundingSphereGenerator.getBoundingSphere(box, sigmaByType, shape);
        walker.setBoundingSphere(bs.center, bs.radius);
        walker.setBoundingSphere(new Vector3D(), 25);
    }

    public void reset() {
        KPlus.E(0);
        KMinus.E(0);
        VPlus.E(0);
        VMinus.E(0);
        hits = 0;
        totalWalks = 0;
    }

    public void setNumWalks(long n) {
        this.numWalks = n;
    }

    protected void recordHit(Vector startPoint, Vector endPoint) {
        Vector boundingSphereCenter = walker.boundingSphereCenter;
        double boundingSphereRadius = walker.boundingSphereRadius;
        char[] walkCharges = new char[3];
        Vector KPlusData = new Vector3D();
        Vector KMinusData = new Vector3D();
        Tensor VPlusData = new Tensor3D();
        Tensor VMinusData = new Tensor3D();
        ++this.hits;

        for(int dim = 0; dim < 3; ++dim) {
            double probability = 0.5 + (startPoint.getX(dim) - boundingSphereCenter.getX(dim)) / (2.0 * boundingSphereRadius);
            if (probability > this.random.nextDouble()) {
                walkCharges[dim] = '+';
                KPlusData.setX(dim, 1.0);
                VPlusData.setComponent(dim, 0, endPoint.getX(0) - boundingSphereCenter.getX(0));
                VPlusData.setComponent(dim, 1, endPoint.getX(1) - boundingSphereCenter.getX(1));
                VPlusData.setComponent(dim, 2, endPoint.getX(2) - boundingSphereCenter.getX(2));
            } else {
                walkCharges[dim] = '-';
                KMinusData.setX(dim, 1.0);
                VMinusData.setComponent(dim, 0, endPoint.getX(0) - boundingSphereCenter.getX(0));
                VMinusData.setComponent(dim, 1, endPoint.getX(1) - boundingSphereCenter.getX(1));
                VMinusData.setComponent(dim, 2, endPoint.getX(2) - boundingSphereCenter.getX(2));
            }
        }

        this.KPlus.PE(KPlusData);
        this.KMinus.PE(KMinusData);
        this.VPlus.PE(VPlusData);
        this.VMinus.PE(VMinusData);
    }

    protected void recordMiss() {
    }

    public void actionPerformed() {
        if (doEq21) {
            reset();
            Vector farthestPoint = BoundingSphereGenerator.findFarthestPoint(box, walker.boundingSphereCenter, sigmaByType, shape);
            double R = walker.boundingSphereRadius;
            double r2 = farthestPoint.Mv1Squared(walker.boundingSphereCenter);
            if (r2 > R*R) {
                throw new RuntimeException("Bounding sphere radius ("+R+") is not large enough, needs to be at least "+Math.sqrt(r2));
            }
        }
        for(long walkNum = 0L; walkNum < this.numWalks; ++walkNum) {
            WalkerExterior.WalkResult result = this.walker.walk();
            if (result.hitObject) {
                this.recordHit(result.startPoint, result.endPoint);
            } else {
                this.recordMiss();
            }
        }

        if (doEq21) {
            double t = hits / (double) this.numWalks;
            Vector u = new Vector3D();
            u.Ev1Mv2(KPlus, KMinus);
            u.TE(1.0 / numWalks);
            Tensor v = new Tensor3D();
            v.E(VPlus);
            v.PE(VMinus);
            v.TE(1.0 / numWalks);
            Tensor w = new Tensor3D();
            w.E(VPlus);
            w.ME(VMinus);
            w.TE(1.0 / numWalks);

            Tensor polarizabililtyTensor = this.computePolarizability(u, v, w, t);
            double q_eta = this.computePadeApproximant(polarizabililtyTensor);
            uAvg.PE(u);
            qvAvg.PEa1Tt1(q_eta, v);
            qwAvg.PEa1Tt1(q_eta, w);

            tAvg += t;
        }


        this.totalWalks += this.numWalks;
        totalSteps++;
    }

    public double[] getAllData() {
        return getAllData(false);
    }

    public double[] getAllData(boolean justW) {
        double t = hits / (double)this.totalWalks;
        double Rh = t * walker.boundingSphereRadius;
//        System.out.println(t+" "+walker.boundingSphereRadius+" "+Rh);

        Tensor polarizabililtyTensor = this.computePolarizability();
        Tensor polarizabililtyTensorJustW = this.computePolarizability(false);
//        System.out.println("polarizability tensor\n" + polarizabililtyTensor);
        double q_eta = doEq21 ? Double.NaN : this.computePadeApproximant(polarizabililtyTensor);
//        System.out.println("q: "+q_eta+"\n");
        // ZENO assumes that error in q_eta is q_eta * 0.015
        double meanPolarizability = polarizabililtyTensor.trace() / 3.0;
        double meanPolarizabilityJustW = polarizabililtyTensorJustW.trace() / 3.0;
        double intrinsicViscosityConventional = q_eta * meanPolarizability / this.box.getLeafList().size();
        double intrinsicViscosityConventionalJustW = q_eta * meanPolarizabilityJustW / this.box.getLeafList().size();
        // ZENO volume computed from interior walker
        double volume = computeVolume();
        double intrinsicConductivity = meanPolarizability / volume;
        double practicalConductivity = meanPolarizability / this.box.getLeafList().size();
        double intrinsicConductivityJustW = meanPolarizabilityJustW / volume;
        double practicalConductivityJustW = meanPolarizabilityJustW / this.box.getLeafList().size();
        // std product uncertainty
        double intrinsicViscosity = q_eta * intrinsicConductivity;
        double intrinsicViscosityJustW = q_eta * intrinsicConductivityJustW;
        if (justW) {
            return new double[]{Rh, meanPolarizabilityJustW, intrinsicConductivityJustW, practicalConductivityJustW, q_eta, intrinsicViscosityJustW, intrinsicViscosityConventionalJustW};
        }
        return new double[]{Rh, meanPolarizability, intrinsicConductivity, practicalConductivity, q_eta, intrinsicViscosity, intrinsicViscosityConventional};
    }

    public double getIntrinsicViscosity() {
        double[] y = getAllData();
        System.out.println("q_eta: " + y[4]);
        System.out.println("mean polarizability: " + y[1]);
        System.out.println("intrinsic conductivity: "+ y[3]);
        return y[5];
    }

    // returns volume occupied by atoms
    public double computeVolume() {
        if (shape == WalkerExterior.Shape.CUBE) return box.getLeafList().size();
        double excludedVolume = 0;
        double nominalVolume = 0;
        for (IAtom a : box.getLeafList()) {
            double rad1 = 0.5*sigmaByType[a.getType().getIndex()];
            nominalVolume += rad1*rad1*rad1;
            for (int j = a.getLeafIndex()+1; j<box.getLeafList().size(); j++) {
                double r2 = a.getPosition().Mv1Squared(box.getLeafList().get(j).getPosition());
                double rad2 = 0.5*sigmaByType[box.getLeafList().get(j).getType().getIndex()];
                double sigma = rad1+rad2;
                double sigma2 = sigma*sigma;
                if (r2 > sigma2) continue;
                double d = Math.sqrt(r2);
                // Google says so
                // https://mathematica.stackexchange.com/questions/73282/how-i-calculate-the-volume-of-multiple-intersecting-spheres
                excludedVolume += Math.PI/(12*d) * Math.pow(sigma - d, 2) * (d*d + 2*d*sigma - 3*Math.pow(rad1-rad2, 2));
            }
        }
        nominalVolume *= 4.0/3.0 * Math.PI;
        return nominalVolume - excludedVolume;
    }

    public double getHydrodynamicRadius() {
        double t = hits / (double)this.totalWalks;
        return t * walker.boundingSphereRadius;
    }

    public Tensor computePolarizability(Vector u, Tensor v, Tensor w, double t) {
        Tensor polarizabilityTensor = new Tensor3D();
//        System.out.println("u: "+u);
//        System.out.println("v: "+v.component(0,0)+" "+v.component(1,1)+" "+v.component(2,2));
//        System.out.println("w\n"+w);

        for(int row = 0; row < 3; ++row) {
            for(int col = 0; col < 3; ++col) {
                double element = 12 * Math.PI * walker.boundingSphereRadius * walker.boundingSphereRadius *
                        (w.component(row, col) - u.getX(row) * v.component(row, col) / t);
                polarizabilityTensor.setComponent(row, col, element);
            }
        }

        Tensor pt = new Tensor3D();
        pt.E(polarizabilityTensor);
        pt.transpose();
        polarizabilityTensor.PE(pt);
        polarizabilityTensor.TE(0.5);
        double l = 1.0; // length scale number
        polarizabilityTensor.TE(Math.pow(l, 3.0));
        return polarizabilityTensor;
    }

    public Tensor computePolarizability() {
        return computePolarizability(true);
    }

    public Tensor computePolarizability(boolean includeUV) {
        double t = hits / (double)this.totalWalks;
        Vector u = new Vector3D();
        Tensor v = new Tensor3D();
        if (!includeUV) {
            u.Ev1Mv2(KPlus, KMinus);
            u.TE(1.0 / totalWalks);
            v.E(VPlus);
            v.PE(VMinus);
            v.TE(1.0 / totalWalks);
        }
        Tensor w = new Tensor3D();
        w.E(VPlus);
        w.ME(VMinus);
        w.TE(1.0 / totalWalks);
//        System.out.println("u: "+u);
//        System.out.println("v: "+v.component(0,0)+" "+v.component(1,1)+" "+v.component(2,2));
//        System.out.println("w: "+w.component(0,0)+" "+w.component(1,1)+" "+w.component(2,2));
//        System.out.println("t: "+t);
        return computePolarizability(u, v, w, t);
    }

    public double computePadeApproximant(Tensor polarizabilityTensor) {
        double alpha1 = 0, alpha2 = 0, alpha3 = 0;
        if (eigenvaluesForPade) {
            Matrix m = new Matrix(polarizabilityTensor.toArray(), 3);
            EigenvalueDecomposition eigen = new EigenvalueDecomposition(m);
            double[] alphas = eigen.getRealEigenvalues();
            alpha1 = alphas[0];
            alpha2 = alphas[1];
            alpha3 = alphas[2];
        }
        else {
            alpha1 = polarizabilityTensor.component(0, 0);
            alpha2 = polarizabilityTensor.component(1, 1);
            alpha3 = polarizabilityTensor.component(2, 2);
            if (alpha2 < alpha1) {
                double t = alpha1;
                alpha1 = alpha2;
                alpha2 = t;
            }

            if (alpha3 < alpha1) {
                double t = alpha1;
                alpha1 = alpha3;
                alpha3 = t;
            }

            if (alpha3 < alpha2) {
                double t = alpha2;
                alpha2 = alpha3;
                alpha3 = t;
            }
        }

        double alpha2_alpha1 = alpha2 / alpha1;
        double alpha3_alpha2 = alpha3 / alpha2;
        if (alpha2_alpha1 < 0.0 || alpha3_alpha2 < 0.0) {
            System.err.println("*** Warning ***\n");
            System.err.println("Could not compute Prefactor for computing intrinsic viscosity, Viscometric radius, or Intrinsic viscosity.  This is likely due to Electric polarizability tensor standard deviations being too high.  Try increasing the number of walks?\n\n");
        }

        double x1 = Math.log(alpha2_alpha1);
        double x2 = Math.log(alpha3_alpha2);
        return this.computePadeApproximant(x1, x2);
    }

    public double computePadeApproximant(double x1, double x2) {
        double[] delta_i = new double[]{4.8, 0.66, -1.247, 0.787};
        double[] k_i = new double[]{0.0, 1.04, 2.012, 2.315};
        double[] b_i = new double[]{0.68, -7.399, 1.048, 0.136};
        double[] t_i = new double[]{0.0, 1.063, 0.895, 4.993};
        double[] B_i = new double[]{1.925, -8.611, 1.652, -0.12};
        double[] q_i = new double[]{0.0, 1.344, 2.029, 1.075};
        double[] c_i = new double[]{13.43, 16.17, 0.51, -5.86};
        double[] r_i = new double[]{0.0, 0.489, 0.879, 2.447};
        double[] A_i = new double[]{16.23, -15.92, 14.83, -3.74};
        double[] v_i = new double[]{0.0, 0.462, 1.989, 4.6};
        double[] m_i = new double[]{2.786, 0.293, -0.11, 0.012};
        double[] u_i = new double[]{0.0, 0.556, 2.034, 3.024};
        double delta = 0.0;
        double b = 0.0;
        double B = 0.0;
        double c = 0.0;
        double A = 0.0;
        double m = 0.0;

        for(int i = 0; i < 4; ++i) {
            delta += delta_i[i] * Math.exp(-k_i[i] * x1);
            b += b_i[i] * Math.exp(-t_i[i] * x1);
            B += B_i[i] * Math.exp(-q_i[i] * x1);
            c += c_i[i] * Math.exp(-r_i[i] * x1);
            A += A_i[i] * Math.exp(-v_i[i] * x1);
            m += m_i[i] * Math.exp(-u_i[i] * x1);
        }

        return (delta * A + c * x2 + b * Math.pow(x2, 2.0) + 4.0 * Math.pow(x2, m)) / (6.0 * A + 6.0 * c * x2 / delta + B * Math.pow(x2, 2.0) + 5.0 * Math.pow(x2, m));
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    double maxR;

    public double getIntrinsicViscosityEq21() {
        Vector uAvgAvg = new Vector3D();
        uAvgAvg.Ea1Tv1(1.0/totalSteps, uAvg);
        Vector qvAvgAvg = new Vector3D();
        for (int i=0; i<qvAvgAvg.getD(); i++) qvAvgAvg.setX(i, qvAvg.component(i, i)/totalSteps);
        Tensor qwAvgAvg = new Tensor3D();
        qwAvgAvg.E(qwAvg);
        qwAvgAvg.TE(1.0/totalSteps);
        double tAvgAvg = tAvg / totalSteps;
        System.out.println("u: "+uAvgAvg);
        System.out.println("qv: "+qvAvgAvg);
        System.out.println("qw: "+qwAvgAvg.component(0,0)+" "+qwAvgAvg.component(1,1)+" "+qwAvgAvg.component(2,2));

        double R = walker.boundingSphereRadius;
        return 4*Math.PI*R*R/box.getLeafList().size() * (qwAvgAvg.trace() - uAvgAvg.dot(qvAvgAvg)/tAvgAvg);
    }

    @Override
    public IData getData() {
        reset();
        resetBoundingSphere();
        double bsr = walker.boundingSphereRadius;
        if (bsr > maxR) {
            maxR = bsr;
            System.out.println("=> "+maxR);
        }
        actionPerformed();
        double[] y = data.getData();
        double[] d = getAllData();
        double[] dJustW = getAllData(true);
        // Rh, meanPolarizability, intrinsicConductivity, practicalConductivity, q_eta, intrinsicViscosity, intrinsicViscosityConventional
        // => hydrodynamic radius, intrinsic conductivity, practical conductivity, intrinsic viscosity, practical viscosity
        y[0] = d[0];
        y[1] = d[2];
        y[2] = d[3];
        y[3] = d[4];
        y[4] = d[5];
        y[5] = d[6];
        y[6] = dJustW[6];
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }


    public MeterHydrodynamicRadius makeMeterHydrodynamicRadius() {
        return new MeterHydrodynamicRadius();
    }

    public class MeterHydrodynamicRadius extends DataSourceScalar {
        public MeterHydrodynamicRadius() {
            super("hydrodynamic radius", Length.DIMENSION);
        }

        @Override
        public double getDataAsScalar() {
            MeterZENO.this.reset();
            MeterZENO.this.resetBoundingSphere();
            MeterZENO.this.actionPerformed();
            return MeterZENO.this.getHydrodynamicRadius();
        }
    }
}


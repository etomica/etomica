/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

/**
 * Computes various quantities that ZENO computes
 */
public class MeterZENO implements IAction {
    protected final Box box;
    protected double shellThickness;
    protected final WalkerExterior walker;
    protected long numWalks;
    protected Vector boundingSphereCenter;
    protected double boundingSphereRadius;
    protected IRandom random;
    protected Vector KPlus;
    protected Vector KMinus;
    protected Tensor VPlus;
    protected Tensor VMinus;
    protected long hits;
    protected long totalWalks;
    protected final double atomRadius = 0.5;

    public MeterZENO(Box box, IRandom random, double[] sigmaByType) {
        this.box = box;
        this.random = random;
        BoundingSphereGenerator.BoundingSphere bs  = BoundingSphereGenerator.getBoundingSphere(box);
        boundingSphereCenter = bs.center;
        this.boundingSphereRadius = bs.radius;
        this.shellThickness = 0.000001 * bs.radius;

        this.walker = new WalkerExterior(box, sigmaByType, random, boundingSphereRadius, boundingSphereCenter, this.shellThickness);
        this.KPlus = new Vector3D();
        this.KMinus = new Vector3D();
        this.VPlus = new Tensor3D();
        this.VMinus = new Tensor3D();
    }

    public void setNumWalks(long n) {
        this.numWalks = n;
    }

    protected void recordHit(Vector startPoint, Vector endPoint) {
        char[] walkCharges = new char[3];
        Vector KPlusData = new Vector3D();
        Vector KMinusData = new Vector3D();
        Tensor VPlusData = new Tensor3D();
        Tensor VMinusData = new Tensor3D();
        ++this.hits;

        for(int dim = 0; dim < 3; ++dim) {
            double probability = 0.5 + (startPoint.getX(dim) - boundingSphereCenter.getX(dim)) / (2.0 * this.boundingSphereRadius);
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
        for(long walkNum = 0L; walkNum < this.numWalks; ++walkNum) {
            WalkerExterior.WalkResult result = this.walker.walk();
            if (result.hitObject) {
                this.recordHit(result.startPoint, result.endPoint);
            } else {
                this.recordMiss();
            }
        }

        this.totalWalks += this.numWalks;
    }

    public double getIntrinsicViscosity() {
        Tensor polarizabililtyTensor = this.computePolarizability();
        System.out.println("polarizability tensor\n" + polarizabililtyTensor);
        double q_eta = this.computePadeApproximant(polarizabililtyTensor);
        // ZENO assumes that error in q_eta is q_eta * 0.015
        System.out.println("q_eta: " + q_eta);
        double meanPolarizability = polarizabililtyTensor.trace() / 3.0;
        System.out.println("mean polarizability: " + meanPolarizability);
        double intrinsicViscosityConventional = q_eta * meanPolarizability / this.box.getLeafList().size();
        // ZENO volume computed from interior walker
        double volume = computeVolume();
        double intrinsicConductivity = meanPolarizability / volume;
        // std product uncertainty
        double intrinsicViscosity = q_eta * intrinsicConductivity;
        return intrinsicViscosity;
    }

    // returns volume occupied by atoms
    public double computeVolume() {
        double excludedVolume = 0;
        for (IAtom a : box.getLeafList()) {
            double rad1 = atomRadius;
            for (int j = a.getLeafIndex()+1; j<box.getLeafList().size(); j++) {
                double r2 = a.getPosition().Mv1Squared(box.getLeafList().get(j).getPosition());
                double rad2 = atomRadius;
                double sigma = rad1+rad2;
                double sigma2 = sigma*sigma;
                if (r2 > sigma2) continue;
                double d = Math.sqrt(r2);
                // Google says so
                // https://mathematica.stackexchange.com/questions/73282/how-i-calculate-the-volume-of-multiple-intersecting-spheres
                excludedVolume += Math.PI/(12*d) * Math.pow(sigma - d, 2) * (d*d + 2*d*sigma - 3*Math.pow(rad1-rad2, 2));
            }
        }
        double nominalVolume = 4.0/3.0 * Math.PI * Math.pow(atomRadius, 3) * box.getLeafList().size();
        return nominalVolume - excludedVolume;
    }

    public double getHydrodynamicRadius() {
        double t = hits / (double)this.totalWalks;
        return t * this.boundingSphereRadius;
    }

    public Tensor computePolarizability() {
        double t = hits / (double)this.totalWalks;
        Vector u = new Vector3D();
        u.Ev1Mv2(KPlus, KMinus);
        u.TE(1.0 / totalWalks);
        Tensor v = new Tensor3D();
        v.E(VPlus);
        v.PE(VMinus);
        v.TE(1.0 / totalWalks);
        Tensor w = new Tensor3D();
        w.E(VPlus);
        w.ME(VMinus);
        w.TE(1.0 / totalWalks);
        Tensor polarizabilityTensor = new Tensor3D();

        for(int row = 0; row < 3; ++row) {
            for(int col = 0; col < 3; ++col) {
                double element = 12 * Math.PI * boundingSphereRadius * boundingSphereRadius *
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

    public double computePadeApproximant(Tensor polarizabilityTensor) {
        double alpha1 = polarizabilityTensor.component(0, 0);
        double alpha2 = polarizabilityTensor.component(1, 1);
        double alpha3 = polarizabilityTensor.component(2, 2);
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
}


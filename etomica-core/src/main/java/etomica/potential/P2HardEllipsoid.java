/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.atom.IAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.spaceNd.VectorND;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

import java.util.List;

/**
 * Evaluates overlap of two identical ellipsoids of revolution at arbitrary separation and orientations.
 * Handles oblate or prolate cases, having a single symmetry axis.
 */
public class P2HardEllipsoid implements IPotential2 {

    protected final Space space;
    protected Vector3D[] a = new Vector3D[4];//these are numbered 1,2,3 because that's how I did it
    protected Vector3D[] b = new Vector3D[4];

    protected final double R1, R2, R3;//half-axes in x, y, z, directions
    protected final double al, be, ga, de, sigmaLarge2, sigmaSmall2;
    protected final double[] A = new double[4];
    protected final double[] B = new double[4];
    protected final Tensor3D Amat = new Tensor3D();
    protected final Tensor3D Bmat = new Tensor3D();
    protected final Tensor3D diag = new Tensor3D();
    protected final Tensor3D Gamma = new Tensor3D();
    protected final Tensor3D C = new Tensor3D();
    protected double w0b, w0c, w1b;
    protected final double w0a, w0d, w1a, w1c, w2a, w2b;
    protected final Vector3D wVec = new Vector3D();//work vector

    /**
     * Constructor with Rperp = 1.0 as default
     * @param space
     * @param Rsym half-length of symmetry axis
     */
    public P2HardEllipsoid(Space space, double Rsym) {
        this(space, Rsym, 1.0);
    }

    /**
     *
     * @param space
     * @param Rsym half-length of symmetry axis
     * @param Rperp half-length of axes perpendicular to symmetry axis
     */
    public P2HardEllipsoid(Space space, double Rsym, double Rperp) {
        this.space = space;
        R1 = Rsym;
        R2 = Rperp;
        R3 = Rperp;
        al = 1 / (R3 * R3);
        ga = 1 / (R1 * R1) - al;
        be = 1 / (R3 * R3);
        de = 1 / (R1 * R1) - be;
        A[1] = B[1] = R1 * R1; //A and B are hardcoded to be the same, but this can be revisited if needed for two ellipsoids of different shapes
        A[2] = B[2] = R2 * R2;
        A[3] = B[3] = R3 * R3;

        w0a = B[1] * B[2] * B[3];
        w0d = A[1] * A[2] * A[3];
        w1a = B[1] * B[2] + B[1] * B[3] + B[2] * B[3];
        w1c = A[1] * A[2] + A[1] * A[3] + A[2] * A[3];
        w2a = B[1] + B[2] + B[3];
        w2b = A[1] + A[2] + A[3];

        double x = 2. * Math.max(R1, R2);
        sigmaLarge2 = x * x; //squared diameter of circumscribed sphere
        x = 2. * Math.min(R1, R2);
        sigmaSmall2 = x * x; //squared diameter of inscribed sphere

        if(Rsym/Rperp < 0.01 || Rsym/Rperp > 200) {
            System.out.println("Overlap may not be reliable for extreme aspect ratios (prolate or oblate)");
            System.out.println("Vieillard-Baron is subject to roundoff errors in this case");
            System.out.println("If application here is needed, consider implementing a proper maximization algorithm with Perram-Wertheim");
            System.out.println("Aspect ratio: "+Rsym/Rperp);
            throw new RuntimeException();
        }
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    /*
    Perram-Wertheim overlap test function.  If this is greater than 1 for any lambda between 0 and 1
    then there is NOT an overlap; else there is.
    Can call with a single value (e.g., lambda = 0.5) to do a prescreen. Otherwise, it may be
    called multiple times as part of a maximization with respect to lambda.
    It is necessary to call setupPW once before calling this function for a new configuration.
     */
    public double FPW(double lam, Vector dr12) {
        double mu = 1 - lam;
        double w0 = lam * lam * lam * w0a + lam * lam * mu * w0b + lam * mu * mu * w0c + mu * mu * mu * w0d;
        double w1 = lam * lam * w1a + lam * mu * w1b + mu * mu * w1c;
        double w2 = lam * w2a + mu * w2b;
        Gamma.E(0.0);
        Gamma.PEa1Tt1(lam, Amat);
        Gamma.PEa1Tt1(mu, Bmat);
        C.E(Gamma);
        C.TE(Gamma);
        C.PEa1Tt1(-w2, Gamma);
        diag.diagE(w1);
        C.PE(diag);
        C.TE(1. / w0);

        wVec.E(dr12);
        C.transform(wVec);
        return lam * mu * dr12.dot(wVec); // lam mu r.C.r
    }

    /*
    Sets up orientation-dependent parameters in Perram-Wertheim method.
    Must be called for new configuration before calling FPW
     */
    public void setupPW(Vector3D a1, Vector3D b1) {
        makeOrthogonalSet(a1, a);
        makeOrthogonalSet(b1, b);

        Amat.Ev1v2(a[1], a[1]);
        Amat.TE(A[1] - A[3]);
        diag.diagE(A[3]);
        Amat.PE(diag);

        Bmat.Ev1v2(b[1], b[1]);
        Bmat.TE(B[1] - B[3]);
        diag.diagE(B[3]);
        Bmat.PE(diag);

        w0b = 0.0;
        w0c = 0.0;
        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
                for (int k = 1; k <= 3; k++) {
                    wVec.E(b[i]);
                    wVec.XE(b[j]);
                    double dot = wVec.dot(a[k]);
                    w0b += B[i] * B[j] * A[k] * dot * dot;

                    wVec.E(a[i]);
                    wVec.XE(a[j]);
                    dot = wVec.dot(b[k]);
                    w0c += A[i] * A[j] * B[k] * dot * dot;

                }
            }
        }
        w0b *= 0.5;
        w0c *= 0.5;

        w1b = w2a * w2b;
        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
                double dot = a[i].dot(b[j]);
                w1b -= A[i] * B[j] * dot * dot;
            }
        }

    }

    // Used by Perram-Wertheim to make arbitrary axis-orientation vectors that are orthogonal to axis of symmetry
    private void makeOrthogonalSet(Vector3D a1, Vector3D[] a) {
        a[1] = new Vector3D(a1);
        a[2] = new Vector3D(-a1.getX(1), a1.getX(0), 0.0);
        a[2].normalize();
        a[3] = new Vector3D(a1.getX(0)*a1.getX(2), a1.getX(1)*a1.getX(2), -(1 - a1.getX(2) * a1.getX(2)) );
        a[3].normalize();
    }

    /*
    Very crude implementation of PW method, using a brute-force sweep of lambda values to identify whether F exceed 1 anywhere
    This is used just to test the PW and VB implementations
    An improved implementation would use an optimization algorithm
     */
    private boolean overlapPW(Vector3D dr, Vector3D a1, Vector3D b1) {
        setupPW(a1, b1);
        double dl = 0.001;
        for(double lam = 0; lam <=1.0; lam+=dl) {
            if (FPW(lam, dr) > 1.0) return false;
        }
        return true;
    }

    /*
    Vieillard-Baron algorithm
     */
    // requires prior call to setupPW
    public boolean overlapVB(Vector3D dr, Vector3D a1, Vector3D b1) {
        double rsq = dr.squared();
        //requires prior call to setupPW for calculation of Amat, Bmat
        double DA = 1/Amat.determinant();
        double DB = 1/Bmat.determinant();
        double aDotb2 = a1.dot(b1);
        aDotb2 *= aDotb2;
        double rDota2 = dr.dot(a1);
        rDota2 *= rDota2;
        double rDotb2 = dr.dot(b1);
        rDotb2 *= rDotb2;


        double PA = -2 * al * be * de - be * ga * de + be * ga * de * aDotb2 -
                3 * al * be*be - be*be*ga + al * DB * rsq - DB + ga * DB * rDota2;
        double PB = -2 * al * be * ga - al * ga * de + al * ga * de * aDotb2 -
                3*al*al*be - al*al*de + be * DA * rsq - DA + de * DA * rDotb2;
        double PAB = al * be * ga * de * (rsq - rsq * aDotb2 + 2*dr.dot(a1)*dr.dot(b1)*a1.dot(b1))
                + al * be*be * ga * rsq + al * be*be * ga * rDota2
                - al * DB * rsq + al*al* be* de * rsq + al*al* be * de * rDotb2 - be * DA * rsq + DA + DB
                + 2*al*al*be*be*rsq - ga*DB*rDota2 - de * DA* rDotb2;

        double n3 = -PA/DB;
        double n2 = -(PA + PB + PAB)/DB;
        double n1 = -PB/DB;

        if((n1 >= 0) && (n2 >= 0) && (n3 >= 0)) return true;

        double n0 = DA/DB;

        double c2 = -0.375 * n3*n3 + n2;
        if(c2 >= 0) return true;

        double c1 = 0.125 * n3*n3*n3 - 0.5*n3*n2 + n1;
        double c0 = -(3./256)* n3*n3*n3*n3 + 0.0625*n3*n3*n2 - 0.25 * n3 * n1 + n0;
        if(c2*c2 - 4*c0 <= 0) return true;

        double term1 = (c2*c2 + 12*c0);
        double term2 = (2*c2*c2*c2 - 72*c2*c0 + 27*c1*c1);
        double D27 = 4*term1*term1*term1 - term2*term2;

        return (D27 < 0);
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        AtomOriented atom1q = (AtomOriented) atom1;
        AtomOriented atom2q = (AtomOriented) atom2;

        double r2 = dr12.squared();

        if (r2 > sigmaLarge2) return 0; // no overlap if further than circumscribed sphere
        if (r2 < sigmaSmall2) return Double.POSITIVE_INFINITY; //overlap if closer than inscribed sphere

        Vector3D aDir = (Vector3D)atom1q.getOrientation().getDirection();
        Vector3D bDir = (Vector3D)atom2q.getOrientation().getDirection();

        setupPW(aDir, bDir);
        if (FPW(0.5, dr12) > 1) return 0; //Perram-Wertheim screening can identify some non-overlaps

        return overlapVB((Vector3D)dr12, aDir, bDir) ? Double.POSITIVE_INFINITY : 0.0;

    }

    /*
    Tests implementation of algorithms by generating random configurations of random ellipsoids and
    checking that all methods agree on whether ellipsoids overlap
     */
    public static void main(String[] args) {
        int[] seeds = RandomNumberGeneratorUnix.getRandSeedArray();
        IRandom random = new RandomMersenneTwister(seeds);

        Vector3D dr = new Vector3D();
        Vector3D a = new Vector3D();
        Vector3D b = new Vector3D();
        boolean agree = true;
        int nTot = 100000;
        int nOverlap = 0;
        for(int i=0; i<nTot; i++) {
            double R = 200 * random.nextDouble() + 0.1;
            double Rperp = 2*random.nextDouble() + 1.0;
            P2HardEllipsoid potential = new P2HardEllipsoid(Space3D.getInstance(),R, Rperp);

            dr.setRandomCube(random);
            dr.TE(15);
            a.setRandomCube(random);
            b.setRandomCube(random);
            a.normalize();
            b.normalize();

            AtomOriented atom1 = new AtomOriented(Space3D.getInstance(), null, true);
            AtomOriented atom2 = new AtomOriented(Space3D.getInstance(), null, true);
            atom1.getOrientation().setDirection(a);
            atom2.getOrientation().setDirection(b);

            boolean overlapPW = potential.overlapPW(dr, a, b);
            boolean overlapVB = potential.overlapVB(dr, a, b);
            boolean overlap = potential.u(dr, atom1, atom2) == 0.0 ? false : true;
            if(overlap) nOverlap++;
            agree = (overlapPW == overlapVB) && (overlapPW == overlap);
            if(!agree) {
                System.out.println("Failure to agree");
                System.out.println("Perram-Wertheim:" + overlapPW);
                System.out.println("Vieillard-Baron:" + overlapVB);
                System.out.println("Combination:\t" + overlap);
                System.out.println("dr: "+dr.toString());
                System.out.println("a: "+a.toString());
                System.out.println("b: "+b.toString());
                System.out.println("R: "+R);
                System.out.println("Rperp: "+Rperp);
                System.out.println(R/Rperp);
                break;
            }
        }
        if(agree) {
            System.out.println("Finished without disagreement");
            System.out.println("Fraction overlap: "+(double)nOverlap/nTot);
        }

    }

}


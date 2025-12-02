package etomica.potential.mbnrg;

import Jama.Matrix;
import etomica.atom.IAtomList;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.models.water.PNWaterGCPM;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.potential.PotentialPolarizable;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Electron;
import etomica.units.Kelvin;

import java.util.Arrays;

import static etomica.math.SpecialFunctions.erfc;

public class PNTTMPolarization implements IPotentialMolecular, PotentialPolarizable {

    protected final double sigma;
    protected final double epsilon, gamma;
    protected final double chargeO, chargeC;
    protected final double core; // = 4.41; //4.41 = 2.1^2; value according to Cummings
    protected final Vector rijVector;
    protected final Vector work, shift;
    protected final Tensor Tunit, Tij;
    protected final double sigmaM;
    protected final double sigmaH;
    protected final double sqrtHMsigmas;
    protected final double massC;
    protected final double massO;
    protected final double a;
    protected final double totalMass;
    protected final double sqrtPiHMsigmas;
    protected final double sqrtPiMMsigmas;
    protected final double alphaC;
    protected final double alphaO;
    protected final double[][] pairPolarization;
    protected Boundary boundary;
    protected Matrix[] Eq, A;
    protected Vector comi, comj;
    private double UpolAtkins;

    public PNTTMPolarization(Space space, Boundary boundary) {
        super();
        this.boundary = boundary;
        a = 2.57;
        sigma = 3.69;
        epsilon = Kelvin.UNIT.toSim(110);
        gamma = 12.75;
        chargeO = Electron.UNIT.toSim(-0.353);
        chargeC = Electron.UNIT.toSim(0.706);
        core = 4.41; //4.41 = 2.1^2; value according to Cummings
        sigmaM = 0.610;
        sigmaH = 0.455;
        sqrtHMsigmas = Math.sqrt(2 * (sigmaH * sigmaH + sigmaM * sigmaM));
        massC = Carbon.INSTANCE.getMass();
        massO = Oxygen.INSTANCE.getMass();
        totalMass = massO * 2 + massC;
        sqrtPiHMsigmas = Math.sqrt(Math.PI * (sigmaH * sigmaH + sigmaM * sigmaM));
        sqrtPiMMsigmas = Math.sqrt(Math.PI * (2 * sigmaM * sigmaM));
        alphaC = 1.47167;
        alphaO = 0.76979;

        comi = space.makeVector();
        comj = space.makeVector();
        shift = space.makeVector();

        rijVector = space.makeVector();

        work = space.makeVector();

        Tunit = space.makeTensor();
        Tunit.E(new double[][]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
        Tij = space.makeTensor();

        Eq = new Matrix[0];
        A = new Matrix[0];
        pairPolarization = new double[10][10];
    }


    public double energy(IMoleculeList molecules) {
        double sum = 0;

        double up = getPolarizationEnergy(molecules);
        if (molecules.size() == 2) {
            int idx0 = molecules.get(0).getIndex();
            int idx1 = molecules.get(1).getIndex();
            if (idx0 > idx1) {
                pairPolarization[idx1][idx0] = up;
            } else {
                pairPolarization[idx0][idx1] = up;
            }
        }
        sum += up;

        return sum;
    }

    public double getLastPolarizationEnergy() {
        return UpolAtkins;
    }

    public double getPolarizationEnergy(IMoleculeList molecules) {

        final int molCount = molecules.size();
        if (Eq.length < molCount + 1) {
            Eq = Arrays.copyOf(Eq, molCount + 1);
            A = Arrays.copyOf(A, molCount + 1);
        }
        if (Eq[molCount] == null) {
            Eq[molCount] = new Matrix(3 * molCount, 1);
            A[molCount] = new Matrix(3 * molCount, 3 * molCount);

            for (int i = 0; i < 3 * molCount; i++) {
                A[molCount].set(i, i, 1);
            }
        }
        final Matrix myEq = Eq[molCount];
        final Matrix myA = A[molCount];
        for (int i = 0; i < 3 * molCount; i++) {
            myEq.set(i, 0, 0);
        }

        /*
         * Finding the Electric fields at the center of mass of each molecule i due to molecule j, Eqi
         * kmb, 8/7/06
         */

        for (int i = 0; i < molecules.size(); i++) {
            IAtomList iLeafAtoms = molecules.get(i).getChildList();
            Vector Cir = iLeafAtoms.get(0).getPosition();
            Vector Oi1r = iLeafAtoms.get(1).getPosition();
            Vector Oi2r = iLeafAtoms.get(2).getPosition();

            comi.Ea1Tv1(massO, Oi1r);
            comi.PEa1Tv1(massO, Oi2r);
            comi.PEa1Tv1(massC, Cir);
            comi.TE(1.0 / totalMass);//c.o.m of molecule i

            for (int j = 0; j < molecules.size(); j++) {
                if (i == j) continue;
                IAtomList jLeafAtoms = molecules.get(j).getChildList();
                Vector Cjr = jLeafAtoms.get(0).getPosition();
                Vector Oj1r = jLeafAtoms.get(1).getPosition();
                Vector Oj2r = jLeafAtoms.get(2).getPosition();

                work.Ev1Mv2(Cir, Cjr);
                shift.Ea1Tv1(-1, work);
                boundary.nearestImage(work);
                shift.PE(work);
                final boolean zeroShift = shift.squared() < 0.1;
                double comtoO1, comtoO2;

                if (zeroShift) {

                    comtoO1 = Math.sqrt(comi.Mv1Squared(Oj1r));
                    comtoO2 = Math.sqrt(comi.Mv1Squared(Oj2r));
                } else {
                    shift.PE(comi);
                    comtoO1 = Math.sqrt(Oj1r.Mv1Squared(shift));
                    shift.ME(comi);

                    shift.PE(comi);
                    comtoO2 = Math.sqrt(Oj2r.Mv1Squared(shift));
                    shift.ME(comi);
                }


                // For molecules that are far apart, fac=chargeX/comWtoX^3, but we add up
                // facs for H and M, which mostly cancel each other out, so we lose quite
                // a bit of precision (~2-3 digits).
                double pol16 = Math.pow(alphaC*alphaO, 1.0/6.0);
                double fac = chargeO / (comtoO1 * comtoO1) *
                        (1 - Math.exp(-a * Math.pow(comtoO1/pol16, 4)));
                work.Ev1Mv2(comi, Oj1r);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i * 3 + 0, 0, myEq.get(i * 3 + 0, 0) + work.getX(0));
                myEq.set(i * 3 + 1, 0, myEq.get(i * 3 + 1, 0) + work.getX(1));
                myEq.set(i * 3 + 2, 0, myEq.get(i * 3 + 2, 0) + work.getX(2));

                fac = chargeO / (comtoO2 * comtoO2) *
                        (1 - Math.exp(-a * Math.pow(comtoO2/pol16, 4)));
                work.Ev1Mv2(comi, Oj2r);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i * 3 + 0, 0, myEq.get(i * 3 + 0, 0) + work.getX(0));
                myEq.set(i * 3 + 1, 0, myEq.get(i * 3 + 1, 0) + work.getX(1));
                myEq.set(i * 3 + 2, 0, myEq.get(i * 3 + 2, 0) + work.getX(2));


                if (i < j) {
                    double OOr2;
                    if (zeroShift) {
                        OOr2 = O1r.Mv1Squared(Ojr);
                    } else {
                        shift.PE(O1r);
                        OOr2 = Ojr.Mv1Squared(shift);
                        shift.ME(O1r);
                    }
                    if (OOr2 < core) {
                        UpolAtkins = Double.NaN;
                        return UpolAtkins;
                    }
                    comj.Ea1Tv1(massH, Hj1r);
                    comj.PEa1Tv1(massO, Ojr);
                    comj.PEa1Tv1(massH, Hj2r);
                    comj.TE(1.0 / totalMass);

                    rijVector.Ev1Mv2(comj, comi);

                    double r12 = Math.sqrt(rijVector.squared());


                    double f = (1 - erfc(r12 / (2 * sigmaM))) - (r12 / (sigmaM * Math.sqrt(Math.PI)) + (r12 * r12 * r12) / (6 * Math.sqrt(Math.PI) * sigmaM * sigmaM * sigmaM)) * Math.exp(-r12 * r12 / (4 * sigmaM * sigmaM));

                    double g = (1 - erfc(r12 / (2 * sigmaM))) - (r12 / (sigmaM * Math.sqrt(Math.PI))) * Math.exp(-r12 * r12 / (4 * sigmaM * sigmaM));

                    // Filling the unit matrix I
                    Tij.Ev1v2(rijVector, rijVector);//Each tensor Tij is a 3X3 matrix

                    Tij.TE(3 * f / (r12 * r12));

                    Tij.PEa1Tt1(-g, Tunit);
                    Tij.TE(1 / (r12 * r12 * r12));

                    //Try matrix inversion solution with Jama library

                    Tij.TE(alphaPol);

                    int mOffset = i * 3;
                    int nOffset = j * 3;
                    for (int m = 0; m < 3; m++) {
                        for (int n = 0; n < 3; n++) {
                            myA.set(mOffset + m, nOffset + n, -Tij.component(m, n));
                            myA.set(nOffset + n, mOffset + m, -Tij.component(n, m));
                        }
                    }
                }
            }

        }
        //x here represents P (almost).
        //For x to be P, the A of the Ax=b actually needs an extra factor of
        //alphaPol.  We'll add that bit in when we calculate UpolAtkins.

        Matrix x = myA.solve(myEq);//myA*x=myEq

        if (false) {
            // this is (mathematically) what we want.  But Jama is slow.
            UpolAtkins = -0.5 * (x.transpose().times(myEq)).get(0, 0) * alphaPol;
        } else {
            UpolAtkins = 0;
            for (int i = 0; i < 3 * molCount; i++) {
                UpolAtkins += x.get(i, 0) * myEq.get(i, 0);
            }
            UpolAtkins *= -0.5 * alphaPol;
        }

        // only needed for more complicated Eq8 from Cummings paper
        if (false) {

            // for the sake of clarity (over perf), just multiply x by alphaPol
            // (see comment above about A lacking alphaPol)
            x.timesEquals(alphaPol);
            Matrix Ep = myA.times(x).minus(x);
            Ep.timesEquals(-1 / alphaPol);

            double x2NormF = x.normF();
            double UpolEquation8 = 2 * UpolAtkins - 0.5 * (x.transpose().times(Ep).get(0, 0)) + (0.5 / alphaPol) * (x2NormF * x2NormF);

            if (Math.abs(UpolAtkins - UpolEquation8) > 1.e-6) {
                throw new RuntimeException("oops " + UpolAtkins + " " + UpolEquation8);
            }
        }

        return UpolAtkins;
    }

}


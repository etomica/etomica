/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import Jama.Matrix;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.potential.IPotentialMolecular;
import etomica.potential.PotentialMolecular;
import etomica.potential.PotentialPolarizable;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.Arrays;

/**
 * GCPM Water potential class.  This class assumes assumes no periodic
 * boundaries exist.  The polarization energy is solved for using matrix
 * inversion rather than iteration, so this class may not be suitable for
 * large systems.
 *
 * @author Ken
 */
public class PNWaterGCPM extends PotentialMolecular implements PotentialPolarizable {
    protected final MoleculePair pair;
    protected final double sigma;
    protected final double epsilon, gamma;
    protected final double chargeH, chargeM;
    protected final double core; // = 4.41; //4.41 = 2.1^2; value according to Cummings
    protected final Vector rijVector;
    protected final Vector work, shift;
    protected final Tensor Tunit, Tij;
    protected final double sigmaM;
    protected final double sigmaH;
    protected final double sqrtHMsigmas;
    protected final double massH;
    protected final double massO;
    protected final double totalMass;
    protected final double sqrtPiHMsigmas;
    protected final double sqrtPiMMsigmas;
    protected final double alphaPol;
    protected final double[][] pairPolarization;
    protected Boundary boundary;
    protected Matrix[] Eq, A;
    protected Vector comWi, comWj;
    protected Component component;
    private double UpolAtkins;

    public PNWaterGCPM(Space space) {
        super(Integer.MAX_VALUE, space);
        pair = new MoleculePair();
        sigma = 3.69;
        epsilon = Kelvin.UNIT.toSim(110);
        gamma = 12.75;
        chargeH = Electron.UNIT.toSim(0.6113);
        chargeM = Electron.UNIT.toSim(-1.2226);
        core = 4.41; //4.41 = 2.1^2; value according to Cummings
        sigmaM = 0.610;
        sigmaH = 0.455;
        sqrtHMsigmas = Math.sqrt(2 * (sigmaH * sigmaH + sigmaM * sigmaM));
        massH = Hydrogen.INSTANCE.getMass();
        massO = Oxygen.INSTANCE.getMass();
        totalMass = massH * 2 + massO;
        sqrtPiHMsigmas = Math.sqrt(Math.PI * (sigmaH * sigmaH + sigmaM * sigmaM));
        sqrtPiMMsigmas = Math.sqrt(Math.PI * (2 * sigmaM * sigmaM));
        alphaPol = 1.444;
        component = Component.FULL;

        comWi = space.makeVector();
        comWj = space.makeVector();
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

    public void setComponent(Component comp) {
        component = comp;
    }

    public double energy(IMoleculeList molecules) {
        double sum = 0;

        if (component != Component.INDUCTION) {
            for (int i = 0; i < molecules.getMoleculeCount() - 1; i++) {
                pair.atom0 = molecules.getMolecule(i);
                for (int j = i + 1; j < molecules.getMoleculeCount(); j++) {
                    pair.atom1 = molecules.getMolecule(j);
                    sum += getNonPolarizationEnergy(pair);
                    if (Double.isInfinite(sum)) {
                        return sum;
                    }
                }
            }
        }
        if (component != Component.TWO_BODY) {
            double up = getPolarizationEnergy(molecules);
            if (molecules.getMoleculeCount() == 2) {
                int idx0 = molecules.getMolecule(0).getIndex();
                int idx1 = molecules.getMolecule(1).getIndex();
                if (idx0 > idx1) {
                    pairPolarization[idx1][idx0] = up;
                } else {
                    pairPolarization[idx0][idx1] = up;
                }
            }
            sum += up;
        }
        return sum;
    }

    public PNWaterGCPMCached makeCachedPairPolarization() {
        return new PNWaterGCPMCached();
    }

    /**
     * This returns the pairwise-additive portion of the GCPM potential for a
     * pair of atoms (dispersion + fixed-charge electrostatics)
     */
    public double getNonPolarizationEnergy(IMoleculeList atoms) {

        IAtomList water1Atoms = atoms.getMolecule(0).getChildList();
        IAtomList water2Atoms = atoms.getMolecule(1).getChildList();

        Vector O1r = water1Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        Vector O2r = water2Atoms.getAtom(SpeciesWater4P.indexO).getPosition();

        work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1, work);
        boundary.nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;

        double r2 = work.squared();

        if (r2 <= core) {
            return Double.POSITIVE_INFINITY;
        }

        Vector H11r = water1Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        Vector H12r = water1Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();
        Vector H21r = water2Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        Vector H22r = water2Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();

        Vector M1r = water1Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        Vector M2r = water2Atoms.getAtom(SpeciesWater4P.indexM).getPosition();

        double r = Math.sqrt(r2);
        double rOverSigma = r / sigma;
        double sigma2OverR2 = 1 / (rOverSigma * rOverSigma);
        double sixOverGamma = 6 / gamma;

        double sum = epsilon / (1 - sixOverGamma) * (sixOverGamma * Math.exp(gamma * (1 - rOverSigma)) - sigma2OverR2 * sigma2OverR2 * sigma2OverR2);//exp-6 potential(Udisp)

        if (zeroShift) {
            r2 = H11r.Mv1Squared(H21r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            r2 = H11r.Mv1Squared(H22r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            r2 = H12r.Mv1Squared(H21r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            r2 = H12r.Mv1Squared(H22r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            r2 = M1r.Mv1Squared(H21r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            r2 = M1r.Mv1Squared(H22r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            r2 = M2r.Mv1Squared(H11r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            r2 = M2r.Mv1Squared(H12r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            r2 = M1r.Mv1Squared(M2r);
            sum += chargeM * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaM)));
        } else {
            shift.PE(H11r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            shift.PE(H11r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            shift.PE(H12r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            shift.PE(H12r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeH * chargeH / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaH)));

            shift.PE(M1r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            shift.PE(M1r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            shift.PE(H11r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(H11r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            shift.PE(H12r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(H12r);
            sum += chargeH * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / sqrtHMsigmas));

            shift.PE(M1r);
            r2 = M2r.Mv1Squared(shift);
            shift.ME(M1r);
            sum += chargeM * chargeM / Math.sqrt(r2) * (1 - SpecialFunctions.erfc(Math.sqrt(r2) / (2 * sigmaM)));

        }
        return sum;
    }

    /**
     * This returns the polarizable portion of the GCPM potential for any
     * number of atoms.
     */
    public double getPolarizationEnergy(IMoleculeList atoms) {

        final int atomCount = atoms.getMoleculeCount();
        if (Eq.length < atomCount + 1) {
            Eq = (Matrix[]) Arrays.resizeArray(Eq, atomCount + 1);
            A = (Matrix[]) Arrays.resizeArray(A, atomCount + 1);
        }
        if (Eq[atomCount] == null) {
            Eq[atomCount] = new Matrix(3 * atomCount, 1);
            A[atomCount] = new Matrix(3 * atomCount, 3 * atomCount);

            for (int i = 0; i < 3 * atomCount; i++) {
                A[atomCount].set(i, i, 1);
            }
        }
        final Matrix myEq = Eq[atomCount];
        final Matrix myA = A[atomCount];
        for (int i = 0; i < 3 * atomCount; i++) {
            myEq.set(i, 0, 0);
        }

        /*
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         */

        for (int i = 0; i < atoms.getMoleculeCount(); i++) {
            IAtomList iLeafAtoms = atoms.getMolecule(i).getChildList();
            Vector O1r = iLeafAtoms.getAtom(SpeciesWater4P.indexO).getPosition();
            Vector H11r = iLeafAtoms.getAtom(SpeciesWater4P.indexH1).getPosition();
            Vector H12r = iLeafAtoms.getAtom(SpeciesWater4P.indexH2).getPosition();

            comWi.Ea1Tv1(massH, H11r);
            comWi.PEa1Tv1(massO, O1r);
            comWi.PEa1Tv1(massH, H12r);
            comWi.TE(1.0 / totalMass);//c.o.m of molecule i

            for (int j = 0; j < atoms.getMoleculeCount(); j++) {
                if (i == j) continue;
                IAtomList jLeafAtoms = atoms.getMolecule(j).getChildList();
                Vector Mjr = jLeafAtoms.getAtom(SpeciesWater4P.indexM).getPosition();
                Vector Ojr = jLeafAtoms.getAtom(SpeciesWater4P.indexO).getPosition();
                Vector Hj1r = jLeafAtoms.getAtom(SpeciesWater4P.indexH1).getPosition();
                Vector Hj2r = jLeafAtoms.getAtom(SpeciesWater4P.indexH2).getPosition();

                work.Ev1Mv2(O1r, Ojr);
                shift.Ea1Tv1(-1, work);
                boundary.nearestImage(work);
                shift.PE(work);
                final boolean zeroShift = shift.squared() < 0.1;
                double comWtoH1, comWtoH2, comWtoM;

                if (zeroShift) {

                    comWtoH1 = Math.sqrt(comWi.Mv1Squared(Hj1r));
                    comWtoH2 = Math.sqrt(comWi.Mv1Squared(Hj2r));
                    comWtoM = Math.sqrt(comWi.Mv1Squared(Mjr));
                } else {
                    shift.PE(comWi);
                    comWtoH1 = Math.sqrt(Hj1r.Mv1Squared(shift));
                    shift.ME(comWi);

                    shift.PE(comWi);
                    comWtoH2 = Math.sqrt(Hj2r.Mv1Squared(shift));
                    shift.ME(comWi);

                    shift.PE(comWi);
                    comWtoM = Math.sqrt(Mjr.Mv1Squared(shift));
                    shift.ME(comWi);

                }


                // For molecules that are far apart, fac=chargeX/comWtoX^3, but we add up
                // facs for H and M, which mostly cancel each other out, so we lose quite
                // a bit of precision (~2-3 digits).
                double fac = chargeH / (comWtoH1 * comWtoH1 * comWtoH1) * ((1 - SpecialFunctions.erfc(comWtoH1 / sqrtHMsigmas))
                        - Math.sqrt(2) * comWtoH1 / sqrtPiHMsigmas * Math.exp(-comWtoH1 * comWtoH1 / (2 * (sigmaM * sigmaM + sigmaH * sigmaH))));
                work.Ev1Mv2(comWi, Hj1r);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i * 3 + 0, 0, myEq.get(i * 3 + 0, 0) + work.getX(0));
                myEq.set(i * 3 + 1, 0, myEq.get(i * 3 + 1, 0) + work.getX(1));
                myEq.set(i * 3 + 2, 0, myEq.get(i * 3 + 2, 0) + work.getX(2));

                fac = chargeH / (comWtoH2 * comWtoH2 * comWtoH2) * ((1 - SpecialFunctions.erfc(comWtoH2 / sqrtHMsigmas))
                        - Math.sqrt(2) * comWtoH2 / sqrtPiHMsigmas * Math.exp(-comWtoH2 * comWtoH2 / (2 * (sigmaM * sigmaM + sigmaH * sigmaH))));
                work.Ev1Mv2(comWi, Hj2r);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i * 3 + 0, 0, myEq.get(i * 3 + 0, 0) + work.getX(0));
                myEq.set(i * 3 + 1, 0, myEq.get(i * 3 + 1, 0) + work.getX(1));
                myEq.set(i * 3 + 2, 0, myEq.get(i * 3 + 2, 0) + work.getX(2));

                fac = chargeM / (comWtoM * comWtoM * comWtoM) * ((1 - SpecialFunctions.erfc(comWtoM / (2 * sigmaM)))
                        - Math.sqrt(2) * comWtoM / sqrtPiMMsigmas * Math.exp(-comWtoM * comWtoM / (4 * sigmaM * sigmaM)));
                work.Ev1Mv2(comWi, Mjr);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i * 3 + 0, 0, myEq.get(i * 3 + 0, 0) + work.getX(0));
                myEq.set(i * 3 + 1, 0, myEq.get(i * 3 + 1, 0) + work.getX(1));
                myEq.set(i * 3 + 2, 0, myEq.get(i * 3 + 2, 0) + work.getX(2));

//                if (i==0) {System.out.println("after "+j); myEq.print(20,12);}

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
                    comWj.Ea1Tv1(massH, Hj1r);
                    comWj.PEa1Tv1(massO, Ojr);
                    comWj.PEa1Tv1(massH, Hj2r);
                    comWj.TE(1.0 / totalMass);

                    rijVector.Ev1Mv2(comWj, comWi);

                    double r12 = Math.sqrt(rijVector.squared());


                    double f = (1 - SpecialFunctions.erfc(r12 / (2 * sigmaM))) - (r12 / (sigmaM * Math.sqrt(Math.PI)) + (r12 * r12 * r12) / (6 * Math.sqrt(Math.PI) * sigmaM * sigmaM * sigmaM)) * Math.exp(-r12 * r12 / (4 * sigmaM * sigmaM));

                    double g = (1 - SpecialFunctions.erfc(r12 / (2 * sigmaM))) - (r12 / (sigmaM * Math.sqrt(Math.PI))) * Math.exp(-r12 * r12 / (4 * sigmaM * sigmaM));

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
            for (int i = 0; i < 3 * atomCount; i++) {
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

    public double getLastPolarizationEnergy() {
        return UpolAtkins;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public enum Component {TWO_BODY, INDUCTION, FULL}

    public class PNWaterGCPMCached implements IPotentialMolecular {

        public double energy(IMoleculeList molecules) {
            int idx0 = molecules.getMolecule(0).getIndex();
            int idx1 = molecules.getMolecule(1).getIndex();
            if (idx0 > idx1) {
                return pairPolarization[idx1][idx0];
            }
            return pairPolarization[idx0][idx1];

        }

        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public void setBox(Box box) {
        }

        public int nBody() {
            return 2;
        }
    }
}

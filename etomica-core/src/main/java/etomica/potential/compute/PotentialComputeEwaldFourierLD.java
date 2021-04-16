package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

import static etomica.math.SpecialFunctions.erfc;
import static java.lang.Math.PI;

/**
 * This computes the Fourier contribution to the Hessian for Ewald sums.  The
 * loops are ordered differently (compared to PotentialCopmuteEwaldFourier) so
 * that the total contribution for each atom pair can be computed here and
 * given to the consumer.
 */
public class PotentialComputeEwaldFourierLD extends PotentialComputeEwaldFourier {

    private void setArraySizes(int numKVectors) {
        if (fExp == null || numKVectors > fExp.length) {
            fExp = new double[numKVectors];
            f6Exp = new double[numKVectors];
        }
    }

    public PotentialComputeEwaldFourierLD(SpeciesManager sm, Box box) {
        super(sm, box);
    }

    /**
     * Returns optimal parameters for the give accuracy target and ratio of
     * computational costs.
     *
     * @param s        accuracy target.  larger values yield more accuracy
     * @param tauRatio ratio of computational cost for pair and Fourier.
     *                 pair cost is for a single pair while Fourier cost is
     *                 for a single atom and single wave vector.  If 0 is
     *                 given, then a default is taken (~2).
     * @return EwaldParams object containing optimal alpha, rCut and kCut for
     * this box
     */
    public EwaldParams getOptimalParams(double s, double tauRatio) {
        int numAtoms = box.getLeafList().size();
        double vol = box.getBoundary().volume();
        // default value based on benchmarks for TestLJMD3DEwald
        if (tauRatio == 0) tauRatio = 33.3 / 17.3;
        EwaldParams params = new EwaldParams();
        params.alpha = Math.pow(tauRatio * Math.pow(PI, 3) * numAtoms / (vol * vol), 1.0 / 6.0);
        params.rCut = s / params.alpha;
        params.kCut = 2 * params.alpha * s;
        return params;
    }

    /**
     * Returns Ewald parameters necessary to achieve accuracy s for the given
     * rCut.
     */
    public EwaldParams getOptimalParamsForCutoff(double s, double rCut) {
        EwaldParams params = new EwaldParams();
        params.rCut = rCut;
        params.alpha = s / rCut;
        params.kCut = 2 * params.alpha * s;
        return params;
    }

    @Override
    public Vector[] getForces() {
        return box.getSpace().makeVectorArray(box.getLeafList().size());
    }

    @Override
    public double getLastVirial() {
        return 0;
    }

    @Override
    public double getLastEnergy() {
        return 0;
    }

    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        if (pc == null || !pc.wantsHessian()) {
            throw new RuntimeException("This is not the PotentialCompute for you");
        }
        int numAtoms = box.getLeafList().size();

        double vol = box.getBoundary().volume();

        Vector bsO2pi = box.getSpace().makeVector();
        bsO2pi.E(box.getBoundary().getBoxSize());
        bsO2pi.TE(0.5/PI);
        int kxMax = (int) (bsO2pi.getX(0) * kCut);
        // cube instead of sphere, so conservatively big
        int[] kMax = {
                kxMax,
                (int) (bsO2pi.getX(1) * kCut),
                (int) (bsO2pi.getX(2) * kCut)
        };
        int[] nk = {kMax[0] + 1, 2 * kMax[1] + 1, 2 * kMax[2] + 1};
        int nkTot = ((2*kxMax + 1) * nk[1] * nk[2] - 1) / 2;

        this.setArraySizes(nkTot);
        this.computeKVectors(kxMax);

        double coeff = 4 * PI / vol;
        double coeffB2 = -2.0 * SQRT_PI * PI * alpha6 * alpha6 * alpha6 / (3.0 * vol);

        for (int ik = 0; ik < this.kxyz2.size(); ik++) {
            double kxyz2 = this.kxyz2.getDouble(ik);
            double expthing = Math.exp(-0.25 * kxyz2 / (alpha * alpha));
            // we could skip this as long as box-length, kCut don't change between calls
            fExp[ik] = coeff * expthing / kxyz2;
            if (alpha6 > 0) {
                double kxyz1 = Math.sqrt(kxyz2);
                double h = kxyz1 / (2 * alpha6);
                double h2 = h * h;
                double exph2 = Math.exp(-h2);
                f6Exp[ik] = coeffB2 * h * h2 * (SQRT_PI * erfc(h) + (0.5 / h2 - 1) / h * exph2);
            }
        }

        IAtomList atoms = box.getLeafList();
        for (int iAtom=0; iAtom<numAtoms; iAtom++) {
            int iType = atoms.get(iAtom).getType().getIndex();
            double qi = chargesByType[iType];
            double Bii = B6[iType][iType];
            if (qi==0 && Bii == 0) {
                continue;
            }

            for (int jAtom=0; jAtom<iAtom; jAtom++) {
                if (pc.skipPair(iAtom, jAtom)) continue;
                int jType = atoms.get(jAtom).getType().getIndex();
                double qj = chargesByType[jType];
                if (qj * qi == 0 && B6[iType][jType] == 0) continue;
                // cos(kr) ?= f(eik.real,ejk.real)

                Vector rij = box.getSpace().makeVector();
                rij.Ev1Mv2(atoms.get(jAtom).getPosition(), atoms.get(iAtom).getPosition());
                Tensor hessianij = box.getSpace().makeTensor();

                for (int ik = 0; ik < this.kxyz2.size(); ik++) {
                    int ikx = this.ik.getInt(ik * 3);
                    int iky = this.ik.getInt(ik * 3 + 1);
                    int ikz = this.ik.getInt(ik * 3 + 2);
                    Vector kvec = Vector.of(ikx, iky, ikz);
                    kvec.DE(bsO2pi);
                    double ijPhiFac = qi*qj*fExp[ik];
                    if (alpha6>0 && B6[iType][jType] != 0) {
                        for (int kB=0; kB<=6; kB++) {
                            ijPhiFac += f6Exp[ik]*b6[iType][kB]*b6[jType][6-kB];
                        }
                    }
                    ijPhiFac *= 2;

                    double kijPhiFac = ijPhiFac * Math.cos(rij.dot(kvec));
                    Tensor phi = box.getSpace().makeTensor();
                    phi.Ev1v2(kvec, kvec);
                    phi.TE(kijPhiFac);
                    hessianij.PE(phi);
                }
                pc.pairComputeHessian(iAtom, jAtom, hessianij);
            }

        }

        return 0;
    }

    @Override
    public double computeOneOld(IAtom atom) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeOne(IAtom atom) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        throw new RuntimeException("nope");
    }

    @Override
    public void processAtomU(double fac) {
        throw new RuntimeException("nope");
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {

            }

            @Override
            public void integratorStepStarted(IntegratorEvent e) {

            }

            @Override
            public void integratorStepFinished(IntegratorEvent e) {

            }
        };
    }
}

package etomica.potential.compute;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.Complex;
import etomica.molecule.IMolecule;
import etomica.potential.IPotential2;
import etomica.potential.P2SoftSphericalSum;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.ewald.P2Ewald1FourierCancel;
import etomica.potential.ewald.P2Ewald6FourierCancel;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import static etomica.math.SpecialFunctions.erfc;
import static etomica.math.SpecialFunctions.factorial;
import static java.lang.Math.PI;

public class PotentialComputeEwaldFourier implements PotentialCompute {

    protected final SpeciesManager sm;
    protected static final double SQRT_PI = Math.sqrt(PI);
    protected final Box box;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected double virialTot = Double.NaN;
    protected boolean includeSelfOne = false;
    protected double uTotSelf, uTotFS = Double.NaN, uOneFS;
    protected double uOneOld, uOneNew;
    protected Vector[] forces;
    protected final int[] atomCountByType;
    protected final Space space;

    protected final double[] chargesByType;
    protected final double[][] B6;
    protected final double[][] b6;
    protected final double[] sigma, epsilon;

    protected final Vector kBasis;

    public void setkCut(double kCut) {
        this.kCut = kCut;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setAlpha6(double alpha6) {
        this.alpha6 = alpha6;
    }

    protected double kCut;
    protected double alpha;
    protected double alpha6;
//    protected Complex[] sFacAtom = new Complex[0]; // Complex for each atom
    protected double[] sFacAtom = new double[0];
//    protected Complex[] sFac = new Complex[0]; // Complex for each kVector
    protected double[] sFac = new double[0];
    // 7 arrays of, Complex for each kVector
    protected final double[][] sFacB = new double[7][0];

    // Array for each spacial dimension, then flattened array of (num kVectors in that dimension)*Complex for each atom
    protected final double[][] eik = new double[3][0];

    protected double[] dsFac = new double[0]; // Complex for each kVector
    // 7 arrays of, Complex for each kVector
    protected final double[][] dsFacB = new double[7][0];
    protected double[] fExp; // double for each kVector
    protected double[] f6Exp; // double for each kVector
    protected int nWaveVectors;
    protected final IntArrayList ik = new IntArrayList();
    protected final DoubleArrayList kxyz2 = new DoubleArrayList();
    protected final Vector3D lastBoxSize = new Vector3D();
    protected double lastKCut = Double.NaN;
    protected boolean didOld;
    public long numMC, tMC;
    public long fTime, fNum;
    boolean first = true;


    private void setArraySizes(int numAtoms, int numKVectors, int[] dimKVectors) {
        if (numAtoms * 2 > sFacAtom.length) {
            sFacAtom = new double[numAtoms * 2];

            forces = new Vector[numAtoms];
            Arrays.setAll(forces, i -> space.makeVector());
        }

        for (int i = 0; i < eik.length; i++) {
            if (2 * dimKVectors[i] * numAtoms > eik[i].length) {
                eik[i] = new double[numAtoms * dimKVectors[i] * 2];
            }
        }

        if (numKVectors * 2 > sFac.length) {
            sFac = new double[numKVectors * 2];

            for (int j = 0; j < sFacB.length; j++) {
                sFacB[j] = new double[numKVectors * 2];
            }

            dsFac = new double[numKVectors * 2];

            for (int j = 0; j < dsFacB.length; j++) {
                dsFacB[j] = new double[numKVectors * 2];
            }

            fExp = new double[numKVectors];
            f6Exp = new double[numKVectors];
        }
    }

    public PotentialComputeEwaldFourier(SpeciesManager sm, Box box) {
        this.sm = sm;
        this.box = box;
        this.space = box.getSpace();
        this.duAtom = new DoubleArrayList(16);
        this.uAtomsChanged = new IntArrayList(16);
        this.forces = new Vector[0];

        int numAtomTypes = sm.getAtomTypeCount();
        this.atomCountByType = new int[numAtomTypes];
        for (ISpecies s : sm.getSpeciesList()) {
            int nMols = box.getNMolecules(s);
            for (AtomType type : s.getAtomTypes()) {
                atomCountByType[type.getIndex()] += nMols;
            }
        }
        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]++;
                }
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]--;
                }
            }
        });

        this.chargesByType = new double[numAtomTypes];
        this.B6 = new double[numAtomTypes][numAtomTypes];
        this.b6 = new double[numAtomTypes][7];
        sigma = new double[numAtomTypes];
        epsilon = new double[numAtomTypes];

        this.kBasis = space.makeVector();
    }

    public static class EwaldParams {
        public double alpha;
        public double rCut;
        public double kCut;

        @Override
        public String toString() {
            return new StringJoiner(", ", EwaldParams.class.getSimpleName() + "[", "]").add("alpha=" + alpha).add("rCut=" + rCut).add("kCut=" + kCut).toString();
        }
    }

    /**
     * Returns a PotentialMasterBonding that computes the contributions that
     * cancel the intramolecular contributions from the Fourier part (computed
     * by this class).  Add the PotentialMasterBonding returned here to your
     * PotentialComputeAggregate.
     *
     * This method must not be called until this class is fully configured.
     * All charges, sigmas, epsilon and alphas must be set.
     */
    public PotentialMasterBonding makeIntramolecularCorrection() {
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        populateIntramolecular(pmBonding);
        return pmBonding;
    }

    /**
     * Populates the given PotentialMasterBonding with pair potentials that
     * cancel the (unwanted) intramolecular contributions computed by this class.
     */
    public void populateIntramolecular(PotentialMasterBonding pmBonding) {
        for (ISpecies species : sm.getSpeciesList()) {
            List<Integer>[] typeList = new ArrayList[sm.getAtomTypeCount()];
            for (int i=0; i<typeList.length; i++) {
                typeList[i] = new ArrayList<>();
            }
            List<AtomType> types = species.getAtomTypes();
            for (int i=0; i<types.size(); i++) {
                typeList[types.get(i).getIndex()].add(i);
            }
            for (int i=0; i<typeList.length; i++) {
                if (typeList[i].size() == 0) continue;
                double qi = chargesByType[i];
                for (int j=i; j<typeList.length; j++) {
                    if (typeList[j].size() == 0 || (i==j && typeList[j].size() == 1)) continue;
                    double qj = chargesByType[j];
                    if (qi*qj==0 && B6[i][j] == 0) continue;
                    List<int[]> pairs = new ArrayList<>();
                    for (int ii = 0; ii < typeList[i].size(); ii++) {
                        int iii = typeList[i].get(ii);
                        if (i==j) {
                            for (int jj=ii+1; jj<typeList[i].size(); jj++) {
                                int jjj = typeList[j].get(jj);
                                if (pmBonding.getBondingInfo().skipBondedChildren(species, iii, jjj)) pairs.add(new int[]{iii,jjj});
                            }
                        }
                        else {
                            for (int jj=0; jj<typeList[j].size(); jj++) {
                                int jjj = typeList[j].get(jj);
                                if (pmBonding.getBondingInfo().skipBondedChildren(species, iii, jjj)) pairs.add(new int[]{iii,jjj});
                            }
                        }
                    }
                    if (pairs.size() == 0) continue;
                    IPotential2 pij = null;
                    if (qi*qj != 0) {
                        pij = new P2Ewald1FourierCancel(qi * qj, alpha);
                    }
                    if (B6[i][j] != 0) {
                        P2Ewald6FourierCancel pij6 = new P2Ewald6FourierCancel(sigma[i], epsilon[i], sigma[j], epsilon[j], alpha6);
                        if (pij != null) {
                            pij = new P2SoftSphericalSum(pij, pij6);
                        }
                        else {
                            pij = pij6;
                        }
                    }
                    pmBonding.setBondingPotentialPair(species, pij, pairs);
                }
            }
        }
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

    public void setCharge(AtomType type, double charge) {
        this.chargesByType[type.getIndex()] = charge;
    }

    public void setR6Coefficient(AtomType type, double sigma, double epsilon) {
        int numAtomTypes = this.B6.length;
        int iType = type.getIndex();
        this.sigma[iType] = sigma;
        this.epsilon[iType] = epsilon;
        double sigmak = 1;
        for (int k=0; k<=6; k++) {
            long ck = factorial(6)/(factorial(6-k)*factorial(k));
            b6[iType][k] = 0.25*sigmak*Math.sqrt(ck*epsilon);
            sigmak *= sigma;
        }
        for (int jType=0; jType<numAtomTypes; jType++) {
            B6[iType][jType] = 0;
            for (int k=0; k<=6; k++) {
                B6[iType][jType] += b6[iType][k]*b6[jType][6-k];
            }
            B6[jType][iType] = B6[iType][jType];
        }
    }

    @Override
    public void init() {

    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    @Override
    public double getLastVirial() {
        return virialTot;
    }

    @Override
    public double getLastEnergy() {
        return uTotSelf + uTotFS;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    protected void computeKVectors(int kxMax) {
        Vector bs = box.getBoundary().getBoxSize();
        if (lastBoxSize.equals(bs) && lastKCut == kCut) {
            return;
        }
        lastBoxSize.E(bs);
        lastKCut = kCut;
        this.ik.clear();
        this.kxyz2.clear();
        double kCut2 = kCut * kCut;
        kBasis.E(2 * PI);
        kBasis.DE(bs);

        for (int ikx=0; ikx<=kxMax; ikx++) {
            double kx = ikx * kBasis.getX(0);
            double kx2 = kx * kx;
            double kyCut2 = kCut2 - kx2;
            boolean xpositive = ikx > 0;
            int kyMax = (int) (0.5 * bs.getX(1) * Math.sqrt(kyCut2) / PI);
            for (int iky = -kyMax; iky <= kyMax; iky++) {
                if (!xpositive && iky < 0) continue;
                boolean ypositive = iky > 0;
                double ky = iky * kBasis.getX(1);
                double kxy2 = kx2 + ky * ky;
                int kzMax = (int) (0.5 * bs.getX(2) * Math.sqrt(kCut2 - kxy2) / PI);
                for (int ikz = -kzMax; ikz <= kzMax; ikz++) {
                    if (!xpositive && !ypositive && ikz <= 0) continue;
                    double kz = ikz * kBasis.getX(2);
                    double kxyz2 = kxy2 + kz * kz;
                    ik.add(ikx);
                    ik.add(iky);
                    ik.add(ikz);
                    this.kxyz2.add(kxyz2);
                }
            }
        }
    }

    private void zeroForces() {
        for (Vector force : this.forces) {
            force.E(0);
        }
    }

    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        int numAtoms = box.getLeafList().size();
        virialTot = 0;
        double q2sum = 0;
        double sumBij = 0;
        double sumBii = 0;
        uTotSelf = uTotFS = 0;

        for (int iType = 0; iType < atomCountByType.length; iType++) {
            int iNum = atomCountByType[iType];
            if (iNum == 0) {
                continue;
            }
            double qi = chargesByType[iType];
            q2sum += iNum * qi * qi;
            double Bii = B6[iType][iType];
            if (alpha6 > 0 && Bii != 0) {
                sumBii += iNum * Bii;
                for (int jType = 0; jType < atomCountByType.length; jType++) {
                    sumBij += iNum * atomCountByType[jType] * B6[iType][jType];
                }
            }
        }
        uTotSelf -= alpha / SQRT_PI * q2sum;
        double vol = box.getBoundary().volume();

        if (sumBii > 0) {
            double alpha63 = alpha6 * alpha6 * alpha6;
            uTotSelf -= SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
            virialTot += 3 * SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
            uTotSelf += alpha63 * (alpha63 / 12) * sumBii;
        }
        double uTot = uTotSelf;

        Vector bs = box.getBoundary().getBoxSize();
        int kxMax = (int) (0.5 * bs.getX(0) / PI * kCut);
        // cube instead of sphere, so conservatively big
        int[] kMax = {
                kxMax,
                (int) (0.5 * bs.getX(1) / PI * kCut),
                (int) (0.5 * bs.getX(2) / PI * kCut)
        };
        int[] nk = {kMax[0] + 1, 2 * kMax[1] + 1, 2 * kMax[2] + 1};
        int nkTot = ((2*kxMax + 1) * nk[1] * nk[2] - 1) / 2;

        this.setArraySizes(numAtoms, nkTot, nk);
        this.computeKVectors(kxMax);
        this.zeroForces();

        // We want exp(i dot(k,r)) for every k and every r
        // then sum over atoms, s(k) = sum[exp(i dot(k,r))] and U(k) = s(k) * s*(k)
        //
        // To get there, we first invoke that
        // exp(i dot(k,r)) = exp(i kx rx) exp(i ky ry) exp(i kz rz)
        // the values of ka (a=x,y,z) will be such that ka = 2 l pi / bs[a]
        // so ea(l=2) = ea(l=1) ea(l=1)
        //    ea(l=3) = ea(l=1) ea(l=2)
        // and so on, where ea(l) = exp(i (2 l pi / bs[a]) ra)
        // See Allen & Tildesley for more details
        // https://dx.doi.org/10.1093/oso/9780198803195.001.0001
        // https://github.com/Allen-Tildesley/examples/blob/master/ewald_module.f90
        IAtomList atoms = box.getLeafList();
        for (int a=0; a<3; a++) {
            if (nkTot == 0) break;
            double fac = 2.0* PI / bs.getX(a);
            for (int iAtom=0; iAtom<numAtoms; iAtom++) {
                int idx = iAtom*nk[a];
                if (a>0) idx += kMax[a];
                Complex.ONE.intoArray(eik[a], idx);
                if (nk[a]==1) continue;
                int iType = atoms.get(iAtom).getType().getIndex();
                if (chargesByType[iType] == 0 && B6[iType][iType] == 0) continue;
                Vector ri = atoms.get(iAtom).getPosition();
                new Complex(Math.cos(fac*ri.getX(a)), Math.sin(fac*ri.getX(a))).intoArray(eik[a], idx+1);
                for (int i=2; i<=kMax[a]; i++) {
                    Complex.fromArray(eik[a], idx + 1).times(Complex.fromArray(eik[a], idx + i - 1))
                            .intoArray(eik[a], idx + i);
                }
                if (a==0) continue;
                for (int i=1; i<=kMax[a]; i++) {
                    Complex.fromArray(eik[a], idx + i).conjugate().intoArray(eik[a], idx - i);
                }
            }
        }

        double coeff = 4 * PI / vol;
        double coeffB2 = -2.0 * SQRT_PI * PI * alpha6 * alpha6 * alpha6 / (3.0 * vol);
        double fourierSum = 0;
        double fourierSum6 = 0;
        double virialSum = 0;
        double virialSum6 = 0;

        boolean debug = false;
        double[] sFacOld = new double[0];
        double[][] sFacBold = new double[7][0];
        if (debug) {
            if (sFacOld.length != sFac.length) {
                sFacOld = new double[sFac.length];
                for (int i = 0; i < sFacB.length; i++) {
                    sFacBold[i] = new double[sFacB[i].length];
                }
            }
            System.arraycopy(sFac, 0, sFacOld, 0, sFac.length);
            for (int i = 0; i < sFacB.length; i++) {
                System.arraycopy(sFacB[i], 0, sFacBold[i], 0, sFacB[i].length);
            }
        }
        for (int kB = 0; kB <= 6; kB++) {
            Arrays.fill(sFacB[kB], 0);
        }
        for (int ik = 0; ik < this.kxyz2.size(); ik++) {
            double kxyz2 = this.kxyz2.getDouble(ik);
            int ikx = this.ik.getInt(ik * 3);
            int iky = this.ik.getInt(ik * 3 + 1);
            int ikz = this.ik.getInt(ik * 3 + 2);
            double expthing = Math.exp(-0.25 * kxyz2 / (alpha * alpha));
            Complex.ZERO.intoArray(sFac, ik);
            // we could skip this as long as box-length, kCut don't change between calls
            fExp[ik] = coeff*expthing/kxyz2;
            double hdf6dh = 0;
            if (alpha6 > 0) {
                double kxyz1 = Math.sqrt(kxyz2);
                double h = kxyz1/(2*alpha6);
                double h2 = h*h;
                double exph2 = Math.exp(-h2);
                f6Exp[ik] = coeffB2*h*h2*(SQRT_PI * erfc(h) + (0.5/h2 - 1)/h*exph2);
                hdf6dh = 3*f6Exp[ik] - 1.5*coeffB2*exph2;
            }

            this.handleKVectorSFac(ik, ikx, iky, ikz, nk[0], nk[1], nk[2], kMax[1], kMax[2], pc);

            double x = Complex.fromArray(sFac, ik).times(Complex.fromArray(sFac, ik).conjugate()).real();
            fourierSum += fExp[ik] * x;
            if (alpha>0) {
                double kdfdk = -(2 + kxyz2 / (2 * alpha * alpha)) * fExp[ik];
                double dfqdV = -fExp[ik] / vol - kdfdk / (3 * vol);
                virialSum += 3*vol*dfqdV*x;
            }
            if (alpha6>0) {
                double df6dV = -f6Exp[ik] / vol - hdf6dh / (3 * vol);
                for (int kB = 0; kB <= 3; kB++) {
                    // the real part of y is the same for kB*conj(6-kB) and for conj(kB)*(6-kB)
                    double y = Complex.fromArray(sFacB[kB], ik).times(Complex.fromArray(sFacB[6 - kB], ik).conjugate()).real();
                    if (kB < 3) y *= 2;
                    fourierSum6 += f6Exp[ik] * y;
                    virialSum6 += 3 * vol * df6dV * y;
                }
            }
            if (doForces) {
                this.handleKVectorForces(ik, ikx, iky, ikz);
            }
        }

        fNum += kxyz2.size() * numAtoms;
        this.nWaveVectors = this.kxyz2.size();
        uTotFS = fourierSum + fourierSum6;
        uTot += uTotFS;

        if (debug && !first) {
            boolean broken = false;
            for (int i = 0; i < sFac.length; i++) {
                if (Math.abs(sFac[i] - sFacOld[i]) > 1e-9) {
                    System.out.println("sFac for " + i + " is off " + sFac[i] + " " + sFacOld[i]);
                    broken = true;
                }
                for (int k = 0; k < sFacB.length; k++) {
                    if (Math.abs(sFacB[k][i] - sFacBold[k][i]) > 1e-5) {
                        System.out.println("sFacB for " + k + " " + i + " is off " + sFacB[k][i] + " " + sFacBold[k][i]);
                        broken = true;
                    }
                }
            }
            if (broken) {
                System.out.flush();
                throw new RuntimeException("oops");
            }
        }

        virialTot += virialSum + virialSum6;
        if (doForces && !sm.isPureAtoms()) {
            virialTot += PotentialCompute.computeVirialIntramolecular(forces, box);
        }
        first = false;

        return uTot;
    }

    private void handleKVectorSFac(int ik, int ikx, int iky, int ikz, int nkx, int nky, int nkz, int kMaxY, int kMaxZ, PotentialCallback pc) {
        long t1 = System.nanoTime();
        IAtomList atoms = box.getLeafList();
        int numAtoms = atoms.size();
        for (int iAtom=0; iAtom<numAtoms; iAtom++) {
            int iType = atoms.get(iAtom).getType().getIndex();
            double qi = chargesByType[iType];
            double Bii = B6[iType][iType];
            if (qi==0 && Bii == 0) {
                Complex.ZERO.intoArray(sFacAtom, iAtom);
                continue;
            }
            Complex iContrib = Complex.fromArray(eik[0], iAtom * nkx + ikx)
                    .times(Complex.fromArray(eik[1], iAtom * nky + kMaxY + iky))
                    .times(Complex.fromArray(eik[2], iAtom * nkz + kMaxZ + ikz));
            iContrib.intoArray(sFacAtom, iAtom);
            Complex.fromArray(sFac, ik).plus(iContrib.times(qi)).intoArray(sFac, ik);
            if (alpha6 > 0) {
                for (int kB=0; kB<=6; kB++) {
                    Complex.fromArray(sFacB[kB], ik).plus(iContrib.times(b6[iType][kB]))
                            .intoArray(sFacB[kB], ik);
                }
            }

            if (pc != null && pc.wantsHessian()) {
                double iphiFac = 2*qi*fExp[ik];
                double[] iphiFac6 = new double[7];
                if (alpha6>0) {
                    for (int kB=0; kB<=6; kB++) {
                        iphiFac6[kB] = 2*f6Exp[ik]*b6[iType][kB];
                    }
                }
                Vector bs = box.getBoundary().getBoxSize();
                Vector kvec = Vector.of(ikx, iky, ikz);
                kvec.TE(2.0*PI);
                kvec.DE(bs);

                for (int jAtom=0; jAtom<iAtom; jAtom++) {
                    if (pc.skipPair(iAtom, jAtom)) continue;
                    int jType = atoms.get(jAtom).getType().getIndex();
                    double qj = chargesByType[jType];
                    if (qj * qi == 0 && iphiFac6[0] * B6[jType][jType] == 0) continue;
                    // cos(kr) ?= f(eik.real,ejk.real)
                    double ijPhiFac = qj * iphiFac;
                    if (alpha6 > 0) {
                        for (int kB = 0; kB <= 6; kB++) {
                            ijPhiFac += b6[jType][6 - kB] * iphiFac6[kB];
                        }
                    }
                    ijPhiFac *= (Complex.fromArray(sFacAtom, iAtom).times(Complex.fromArray(sFacAtom, jAtom).conjugate())).real();
                    Tensor phi = box.getSpace().makeTensor();
                    phi.Ev1v2(kvec, kvec);
                    phi.TE(ijPhiFac);
                    pc.pairComputeHessian(iAtom, jAtom, phi);
                }
            }

        }
        fTime += System.nanoTime() - t1;
    }

    private void handleKVectorForces(int ik, int ikx, int iky, int ikz) {
        IAtomList atoms = box.getLeafList();
        int numAtoms = atoms.size();
        double kx = kBasis.getX(0) * ikx;
        double ky = kBasis.getX(1) * iky;
        double kz = kBasis.getX(2) * ikz;
        for (int iAtom=0; iAtom<numAtoms; iAtom++) {
            int iType = atoms.get(iAtom).getType().getIndex();
            double coeffki = alpha==0 ? 0 :
                    2*fExp[ik]*chargesByType[iType] *
                            Complex.fromArray(sFacAtom,iAtom).times(Complex.fromArray(sFac, ik).conjugate()).imaginary();
            double coeffki6 = 0;
            if (alpha6 > 0) {
                for (int kB = 0; kB <= 6; kB++) {
                    coeffki6 += 2 * f6Exp[ik] * b6[iType][kB] *
                            Complex.fromArray(sFacAtom, iAtom).times(Complex.fromArray(sFacB[6 - kB], ik).conjugate()).imaginary();
                }
            }
            coeffki += coeffki6;
            forces[iAtom].PE(Vector.of(coeffki * kx, coeffki * ky, coeffki * kz));
        }
    }

    @Override
    public double computeOneOld(IAtom atom) {
        didOld = true;
        Arrays.fill(dsFac, 0);
        for (int i = 0; i < dsFacB.length; i++) Arrays.fill(dsFacB[i], 0);
        return computeOneInternal(atom, true);
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        didOld = true;
        Arrays.fill(dsFac, 0);
        for (int i = 0; i < dsFacB.length; i++) Arrays.fill(dsFacB[i], 0);
        return computeManyAtomsInternal(true, ((AtomArrayList) molecule.getChildList()).toArray());
    }

    @Override
    public double computeOne(IAtom atom) {
        if (!didOld) return 0;
        didOld = false;
        return computeOneInternal(atom, false);
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        return computeManyAtoms(((AtomArrayList) molecule.getChildList()).toArray());
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        didOld = true;
        Arrays.fill(dsFac, 0);
        for (int i = 0; i < dsFacB.length; i++) Arrays.fill(dsFacB[i], 0);
        return computeManyAtomsInternal(true, atoms);
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        if (!didOld) return 0;
        didOld = false;
        return computeManyAtomsInternal(false, atoms);
    }

    private double computeOneInternal(IAtom atom, boolean oldEnergy) {
        int iAtom = atom.getLeafIndex();
        int iType = atom.getType().getIndex();
        double qi = chargesByType[iType];
        if (oldEnergy) qi = -qi;
        if (qi == 0 && B6[iType][iType] == 0) return 0;

        double uOne = 0;
        if (includeSelfOne) {
            double q2Sum = qi * qi;
            uOne -= alpha / SQRT_PI * q2Sum;
            double vol = box.getBoundary().volume();
            if (alpha6 > 0 && B6[iType][iType] != 0) {
                double sumBii = B6[iType][iType];
                double sumBij = 0;
                for (int jType = 0; jType < atomCountByType.length; jType++) {
                    if (iType != jType) {
                        //  Ni*Nj - (Ni-1)*Nj = Nj
                        sumBij += atomCountByType[jType] * B6[iType][jType];
                    } else {
                        // N*N - (N-1)*(N-1) = 2N-1
                        sumBij += (2 * atomCountByType[jType] - 1) * B6[iType][jType];
                    }
                }
                double alpha63 = alpha6 * alpha6 * alpha6;
                uOne -= SQRT_PI * PI / (6 * vol) * alpha63 * sumBij;
                uOne += (1.0 / 12.0) * alpha63 * alpha63 * sumBii;
            }
        }
        Vector bs = box.getBoundary().getBoxSize();
        int kxMax = (int) (0.5 * bs.getX(0) / PI * kCut);
        // cube instead of sphere, so conservatively big
        int[] kMax = {
                kxMax,
                (int) (0.5 * bs.getX(1) / PI * kCut),
                (int) (0.5 * bs.getX(2) / PI * kCut)
        };
        int[] nk = {kMax[0] + 1, 2 * kMax[1] + 1, 2 * kMax[2] + 1};
        // we could save eik too... would need to update after accepted move, revert after rejection
        Vector ri = atom.getPosition();
        for (int a=0; a<3; a++) {
            double fac = 2.0*PI/bs.getX(a);
            int idx = iAtom*nk[a];
            if (a>0) idx += kMax[a];
            Complex.ONE.intoArray(eik[a], idx);
            new Complex(Math.cos(fac * ri.getX(a)), Math.sin(fac * ri.getX(a)))
                    .intoArray(eik[a], idx + 1);
            for (int i=2; i<=kMax[a]; i++) {
                Complex.fromArray(eik[a], idx + 1).times(Complex.fromArray(eik[a], idx + i - 1))
                        .intoArray(eik[a], idx + i);
            }
            if (a==0) continue;
            for (int i=1; i<=kMax[a]; i++) {
                Complex.fromArray(eik[a], idx + i).conjugate().intoArray(eik[a], idx - i);
            }
        }

        double fourierSum = 0;
        double fourierSum6 = 0;
        for (int ik = 0; ik < this.kxyz2.size(); ik++) {
            int ikx = this.ik.getInt(ik * 3);
            int iky = this.ik.getInt(ik * 3 + 1);
            int ikz = this.ik.getInt(ik * 3 + 2);


            Complex iContrib = Complex.fromArray(eik[0], iAtom * nk[0] + ikx)
                    .times(Complex.fromArray(eik[1], iAtom * nk[1] + kMax[1] + iky))
                    .times(Complex.fromArray(eik[2], iAtom * nk[2] + kMax[2] + ikz));
            if (alpha > 0) {
                Complex.fromArray(dsFac, ik).plus(iContrib.times(qi)).intoArray(dsFac, ik);
            }
            if (alpha6 > 0) {
                for (int kB = 0; kB <= 6; kB++) {
                    Complex.fromArray(dsFacB[kB], ik).plus(iContrib.times((oldEnergy ? -1 : +1) * b6[iType][kB]))
                            .intoArray(dsFacB[kB], ik);
                }
            }

            if (oldEnergy) {
                if (alpha > 0) {
                    Complex sFacMinus = Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik));
                    fourierSum += fExp[ik] * sFacMinus.times(sFacMinus.conjugate()).real();
                }

                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 3; kB++) {
                        Complex sFacMinus = Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacMinus2 = Complex.fromArray(sFacB[6 - kB], ik).plus(Complex.fromArray(dsFacB[6 - kB], ik));
                        fourierSum6 += (kB < 3 ? 2 : 1) * f6Exp[ik] * sFacMinus.times(sFacMinus2.conjugate()).real();
                    }
                }
            } else {
                if (alpha > 0) {
                    Complex sFacNew = Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik));
                    fourierSum += fExp[ik] * sFacNew.times(sFacNew.conjugate()).real();
                }
                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 3; kB++) {
                        Complex sFacNew = Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacNew2 = Complex.fromArray(sFacB[6 - kB], ik).plus(Complex.fromArray(dsFacB[6 - kB], ik));
                        fourierSum6 += (kB < 3 ? 2 : 1) * f6Exp[ik] * sFacNew.times(sFacNew2.conjugate()).real();
                    }
                }
            }
        }

        if (oldEnergy) {
            // uOne = uTot - (energy without this atom)
            uOneFS = uTotFS - fourierSum - fourierSum6;
            uOne += uOneFS;
        } else {
            // (energy with new coordinates) - (energy without the atom)
            uOne += fourierSum + fourierSum6 - (uTotFS - uOneFS);
        }
        return uOne;
    }

    private double computeManyAtomsInternal(boolean oldEnergy, IAtom... atoms) {

        double q2sum = 0;
        double sumBij = 0;
        double sumBii = 0;
        double uOne = 0;
        int[] typesRemoved = new int[atomCountByType.length];

        if (includeSelfOne) {
            for (IAtom atom : atoms) {
                int iType = atom.getType().getIndex();
                double qi = chargesByType[iType];
                q2sum += qi * qi;
                double Bii = B6[iType][iType];
                if (alpha6 > 0 && Bii != 0) {
                    sumBii += Bii;
                    typesRemoved[iType]++;
                }
            }

            if (q2sum == 0 && sumBii == 0) return 0;

            uOne -= alpha / SQRT_PI * q2sum;
            double vol = box.getBoundary().volume();

            if (sumBii > 0) {
                for (int iType = 0; iType < atomCountByType.length; iType++) {
                    int ni = atomCountByType[iType];
                    for (int jType = 0; jType < atomCountByType.length; jType++) {
                        int nj = atomCountByType[jType];
                        sumBij += (ni * nj - (ni - typesRemoved[iType]) * (nj - typesRemoved[jType])) * B6[iType][jType];
                    }
                }

                double alpha63 = alpha6 * alpha6 * alpha6;
                uOne -= SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
                virialTot += 3 * SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
                uOne += alpha63 * (alpha63 / 12) * sumBii;
            }
        }

        Vector bs = box.getBoundary().getBoxSize();
        int kxMax = (int) (0.5 * bs.getX(0) / PI * kCut);
        // cube instead of sphere, so conservatively big
        int[] kMax = {
                kxMax,
                (int) (0.5 * bs.getX(1) / PI * kCut),
                (int) (0.5 * bs.getX(2) / PI * kCut)
        };
        int[] nk = {kMax[0] + 1, 2 * kMax[1] + 1, 2 * kMax[2] + 1};

        // We want exp(i dot(k,r)) for every k and every r
        // then sum over atoms, s(k) = sum[exp(i dot(k,r))] and U(k) = s(k) * s*(k)
        //
        // To get there, we first invoke that
        // exp(i dot(k,r)) = exp(i kx rx) exp(i ky ry) exp(i kz rz)
        // the values of ka (a=x,y,z) will be such that ka = 2 l pi / bs[a]
        // so ea(l=2) = ea(l=1) ea(l=1)
        //    ea(l=3) = ea(l=1) ea(l=2)
        // and so on, where ea(l) = exp(i (2 l pi / bs[a]) ra)
        // See Allen & Tildesley for more details
        // https://dx.doi.org/10.1093/oso/9780198803195.001.0001
        // https://github.com/Allen-Tildesley/examples/blob/master/ewald_module.f90
        for (int a = 0; a < 3; a++) {
            double fac = 2.0 * PI / bs.getX(a);
            for (IAtom atom : atoms) {
                int iAtom = atom.getLeafIndex();
                int idx = iAtom * nk[a];
                if (a > 0) idx += kMax[a];
                Complex.ONE.intoArray(eik[a], idx);
                if (nk[a] == 1) continue;
                int iType = atom.getType().getIndex();
                if (chargesByType[iType] == 0 && B6[iType][iType] == 0) continue;
                Vector ri = atom.getPosition();
                new Complex(Math.cos(fac * ri.getX(a)), Math.sin(fac * ri.getX(a))).intoArray(eik[a], idx + 1);
                for (int i = 2; i <= kMax[a]; i++) {
                    Complex.fromArray(eik[a], idx + 1).times(Complex.fromArray(eik[a], idx + i - 1))
                            .intoArray(eik[a], idx + i);
                }
                if (a == 0) continue;
                for (int i = 1; i <= kMax[a]; i++) {
                    Complex.fromArray(eik[a], idx + i).conjugate().intoArray(eik[a], idx - i);
                }
            }
        }

        numMC += kxyz2.size() * atoms.length;
        long t1 = System.nanoTime();

        double fourierSum = 0;
        double fourierSum6 = 0;
        for (int ik = 0; ik < this.kxyz2.size(); ik++) {
            int ikx = this.ik.getInt(ik * 3);
            int iky = this.ik.getInt(ik * 3 + 1);
            int ikz = this.ik.getInt(ik * 3 + 2);

            for (IAtom atom : atoms) {
                int iAtom = atom.getLeafIndex();
                int iType = atom.getType().getIndex();
                double qi = chargesByType[iType];
                if (oldEnergy) qi = -qi;
                Complex iContrib = Complex.fromArray(eik[0], iAtom * nk[0] + ikx)
                        .times(Complex.fromArray(eik[1], iAtom * nk[1] + kMax[1] + iky))
                        .times(Complex.fromArray(eik[2], iAtom * nk[2] + kMax[2] + ikz));
                if (alpha > 0) {
                    Complex.fromArray(dsFac, ik).plus(iContrib.times(qi)).intoArray(dsFac, ik);
                }
                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 6; kB++) {
                        Complex.fromArray(dsFacB[kB], ik).plus(iContrib.times((oldEnergy ? -1 : +1) * b6[iType][kB]))
                                .intoArray(dsFacB[kB], ik);
                    }
                }
            }

            if (oldEnergy) {
                if (alpha > 0) {
                    Complex sFacMinus = Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik));
                    fourierSum += fExp[ik] * sFacMinus.times(sFacMinus.conjugate()).real();
                }

                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 3; kB++) {
                        Complex sFacMinus = Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacMinus2 = Complex.fromArray(sFacB[6 - kB], ik).plus(Complex.fromArray(dsFacB[6 - kB], ik));
                        fourierSum6 += (kB < 3 ? 2 : 1) * f6Exp[ik] * sFacMinus.times(sFacMinus2.conjugate()).real();
                    }
                }
            } else {
                if (alpha > 0) {
                    Complex sFacNew = Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik));
                    fourierSum += fExp[ik] * sFacNew.times(sFacNew.conjugate()).real();
                }
                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 3; kB++) {
                        Complex sFacNew = Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacNew2 = Complex.fromArray(sFacB[6 - kB], ik).plus(Complex.fromArray(dsFacB[6 - kB], ik));
                        fourierSum6 += (kB < 3 ? 2 : 1) * f6Exp[ik] * sFacNew.times(sFacNew2.conjugate()).real();
                    }
                }
            }
        }
        tMC += System.nanoTime() - t1;
        if (oldEnergy) {
            // uOne = uTot - (energy without this molecule)
            uOneFS = uTotFS - fourierSum - fourierSum6;
            uOne += uOneFS;
            uOneOld = uOne;
        } else {
            // (energy with new coordinates) - (energy without the molecule)
            uOne += fourierSum + fourierSum6 - (uTotFS - uOneFS);
            uOneNew = uOne;
        }
        return uOne;
    }

    @Override
    public void processAtomU(double fac) {
        for (int ik = 0; ik < nWaveVectors; ik++) {
            if (fac == 1) {
                Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik)).intoArray(sFac, ik);
                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 6; kB++) {
                        Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik)).intoArray(sFacB[kB], ik);
                    }
                }
            }
        }
        Arrays.fill(dsFac, 0, nWaveVectors * 2, 0);
        if (alpha6 > 0) {
            for (int kB = 0; kB <= 6; kB++) {
                Arrays.fill(dsFacB[kB], 0, nWaveVectors * 2, 0);
            }
        }
        if (fac == 1) {
            uTotFS += uOneNew - uOneOld;
        }
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

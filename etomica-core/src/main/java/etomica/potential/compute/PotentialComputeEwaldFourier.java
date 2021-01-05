package etomica.potential.compute;

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
import etomica.potential.BondingInfo;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;
import java.util.StringJoiner;

import static etomica.math.SpecialFunctions.erfc;
import static etomica.math.SpecialFunctions.factorial;
import static java.lang.Math.PI;

public class PotentialComputeEwaldFourier implements PotentialCompute {
    private static final double SQRT_PI = Math.sqrt(PI);
    private final Box box;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected double virialTot;
    protected Vector[] forces;
    protected final int[] atomCountByType;
    protected final Space space;
    private final BondingInfo bondingInfo;

    private final double[] chargesByType;
    private final double[][] B6;
    private final double[][] b6;

    private final Vector kBasis;

    public void setkCut(double kCut) {
        this.kCut = kCut;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public void setAlpha6(double alpha6) {
        this.alpha6 = alpha6;
    }

    private double kCut;
    private double alpha;
    private double alpha6;
//    private Complex[] sFacAtom = new Complex[0]; // Complex for each atom
    private double[] sFacAtom = new double[0];
//    private Complex[] sFac = new Complex[0]; // Complex for each kVector
    private double[] sFac = new double[0];
//    private final Complex[][] sFacB = new Complex[7][0]; // 7 arrays of, Complex for each kVector
    private final double[][] sFacB = new double[7][0];

    // Array for each spacial dimension, then flattened array of (num kVectors in that dimension)*Complex for each atom
//    private final Complex[][] eik = new Complex[3][0];
    private final double[][] eik = new double[3][0];

    private double[] dsFac = new double[0]; // Complex for each kVector
//    private final Complex[][] dsFacB = new Complex[7][0]; // 7 arrays of, Complex for each kVector
    private final double[][] dsFacB = new double[7][0];
    private double[] fExp; // double for each kVector
    private double[] f6Exp; // double for each kVector
    private int nWaveVectors;
    private final IntArrayList ik = new IntArrayList();
    private final DoubleArrayList kxyz2 = new DoubleArrayList();
    private final Vector3D lastBoxSize = new Vector3D();
    private double lastKCut = Double.NaN;


    private void setArraySizes(int numAtoms, int numKVectors, int[] dimKVectors) {
        if (numAtoms > sFacAtom.length) {
            sFacAtom = new double[numAtoms*2];

            forces = new Vector[numAtoms];
            Arrays.setAll(forces, i -> space.makeVector());
        }

        for (int i = 0; i < eik.length; i++) {
            if (dimKVectors[i] * numAtoms > eik[i].length) {
                eik[i] = new double[numAtoms * dimKVectors[i] * 2];
            }
        }

        if (numKVectors > sFac.length) {
            sFac = new double[numKVectors*2];

            for (int j = 0; j < sFacB.length; j++) {
                sFacB[j] = new double[numKVectors*2];
            }

            dsFac = new double[numKVectors*2];

            for (int j = 0; j < dsFacB.length; j++) {
                dsFacB[j] = new double[numKVectors*2];
            }

            fExp = new double[numKVectors];
            f6Exp = new double[numKVectors];
        }
    }

    public PotentialComputeEwaldFourier(Simulation sim, Box box, BondingInfo bondingInfo) {
        this.box = box;
        this.space = box.getSpace();
        this.bondingInfo = bondingInfo;
        this.duAtom = new DoubleArrayList(16);
        this.uAtomsChanged = new IntArrayList(16);
        this.forces = new Vector[0];

        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        int numAtomTypes = lastTypeIndex + 1;
        this.atomCountByType = new int[numAtomTypes];
        for (ISpecies s : sim.getSpeciesList()) {
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

        this.kBasis = space.makeVector();
    }

    public static class EwaldParams {
        public double alpha;
        public double rCut;
        public double kCut;

        @Override
        public String toString() {
            return new StringJoiner(", ", EwaldParams.class.getSimpleName() + "[", "]") .add("alpha=" + alpha) .add("rCut=" + rCut) .add("kCut=" + kCut) .toString();
        }
    }

    public EwaldParams getOptimalParams(double s, double tauRatio) {
        int numAtoms = box.getLeafList().size();
        double vol = box.getBoundary().volume();
        // based on simple benchmarks for etomica-cpp
        if (tauRatio == 0) tauRatio = 9;
        EwaldParams params = new EwaldParams();
        params.alpha = Math.pow(tauRatio * Math.pow(PI, 3) * numAtoms / (vol*vol), 1.0/6.0);
        params.rCut = s/params.alpha;
        params.kCut = 2 * params.alpha*s;
        return params;
    }

    public void setCharge(AtomType type, double charge) {
        this.chargesByType[type.getIndex()] = charge;
    }

    public void setR6Coefficient(AtomType type, double sigma, double epsilon) {
        int numAtomTypes = this.B6.length;
        int iType = type.getIndex();
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
    public double getOldEnergy() {
        double q2sum = 0;
        double sumBij = 0;
        double sumBii = 0;
        double uTot = 0;
        for (int iType = 0; iType < atomCountByType.length; iType++) {
            int iNum = atomCountByType[iType];
            if (iNum == 0) { continue; }
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
        uTot -= alpha / SQRT_PI * q2sum;
        double vol = box.getBoundary().volume();

        if (sumBii > 0) {
            double alpha63 = alpha6 * alpha6 * alpha6;
            uTot -= SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
            uTot += alpha63 * (alpha63 / 12) * sumBii;
        }

        for (IMolecule molecule : box.getMoleculeList()) {
            uTot += computeFourierIntramolecular(molecule, false);
        }

        double fourierSum = 0, fourierSum6 = 0;
        for (int ik = 0; ik < nWaveVectors; ik++) {
            double x = Complex.fromArray(sFac, ik).times(Complex.fromArray(sFac, ik).conjugate()).real();
            fourierSum += fExp[ik] * x;
            if (alpha6 > 0) {
                for (int kB = 0; kB <= 6; kB++) {
                    double y = Complex.fromArray(sFacB[kB], ik).times(Complex.fromArray(sFacB[6 - kB], ik).conjugate()).real();
                    fourierSum6 += f6Exp[ik] * y;
                }
            }
        }
        uTot += fourierSum;
        uTot += fourierSum6;
        return uTot;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    private void computeKVectors(int kxMax) {
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
    public double computeAll(boolean doForces) {
        int numAtoms = box.getLeafList().size();
        virialTot = 0;
        double q2sum = 0;
        double sumBij = 0;
        double sumBii = 0;
        double uTot = 0;

        for (int iType = 0; iType < atomCountByType.length; iType++) {
            int iNum = atomCountByType[iType];
            if (iNum == 0) { continue; }
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
        uTot -= alpha / SQRT_PI * q2sum;
        double vol = box.getBoundary().volume();

        if (sumBii > 0) {
            double alpha63 = alpha6 * alpha6 * alpha6;
            uTot -= SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
            virialTot += 3 * SQRT_PI * PI * alpha63 / (6 * vol) * sumBij;
            uTot += alpha63 * (alpha63 / 12) * sumBii;
        }

        for (IMolecule molecule : box.getMoleculeList()) {
            uTot += computeFourierIntramolecular(molecule, doForces);
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
        double coeffB2 = -2.0* SQRT_PI * PI * alpha6*alpha6*alpha6/(3.0*vol);
        double fourierSum = 0;
        double fourierSum6 = 0;
        double virialSum = 0;
        double virialSum6 = 0;

        for (int kB=0; kB<=6; kB++) {
            Arrays.fill(sFacB[kB], 0);
        }
        for (int ik = 0; ik < this.kxyz2.size(); ik++) {
            double kxyz2 = this.kxyz2.getDouble(ik);
            int ikx = this.ik.getInt(ik * 3);
            int iky = this.ik.getInt(ik * 3 + 1);
            int ikz = this.ik.getInt(ik * 3 + 2);
            double expthing = Math.exp(-0.25*kxyz2/(alpha*alpha));
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

            this.handleKVectorSFac(ik, ikx, iky, ikz, nk[0], nk[1], nk[2], kMax[1], kMax[2]);

            double x = Complex.fromArray(sFac, ik).times(Complex.fromArray(sFac, ik).conjugate()).real();
            fourierSum += fExp[ik] * x;
            if (alpha>0) {
                double kdfdk = -(2 + kxyz2 / (2 * alpha * alpha)) * fExp[ik];
                double dfqdV = -fExp[ik] / vol - kdfdk / (3 * vol);
                virialSum += 3*vol*dfqdV*x;
            }
            if (alpha6>0) {
                double df6dV = -f6Exp[ik] / vol - hdf6dh / (3 * vol);
                for (int kB=0; kB<=6; kB++) {
                    double y = Complex.fromArray(sFacB[kB], ik).times(Complex.fromArray(sFacB[6 - kB], ik).conjugate()).real();
                    fourierSum6 += f6Exp[ik] * y;
                    virialSum6 += 3*vol*df6dV*y;
                }
            }
            if (doForces) {
                this.handleKVectorForces(ik, ikx, iky, ikz);
            }
        }

        this.nWaveVectors = this.kxyz2.size();
        uTot += fourierSum + fourierSum6;
        virialTot += virialSum + virialSum6;

        return uTot;
    }

    private void handleKVectorSFac(int ik, int ikx, int iky, int ikz, int nkx, int nky, int nkz, int kMaxY, int kMaxZ) {
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
        }
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

    protected double computeFourierIntramolecular(IMolecule molecule, boolean doForces) {
        IAtomList atoms = molecule.getChildList();
        if (atoms.size() == 1) return 0;
        double twoosqrtpi = 2.0 / SQRT_PI;
        double uTot = 0;
        for (int iAtom = 0; iAtom < atoms.size() - 1; iAtom++) {
            IAtom atom = atoms.get(iAtom);
            int iType = atom.getType().getIndex();
            double qi = chargesByType[iType];
            double Bii = B6[iType][iType];
            if (qi == 0 && Bii == 0) continue;
            Vector ri = atom.getPosition();
            for (int jAtom = iAtom + 1; jAtom < atoms.size(); jAtom++) {
                int jType = atoms.get(jAtom).getType().getIndex();
                double qj = chargesByType[jType];
                double Bij = B6[iType][jType];
                double qiqj = qi * qj;
                if (qiqj == 0 && Bij == 0) continue;
                Vector rj = atoms.get(jAtom).getPosition();
                Vector dr = space.makeVector();
                dr.Ev1Mv2(rj, ri);
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                double r = Math.sqrt(r2);
                double ec = 0;
                if (alpha > 0) {
                    ec = erfc(alpha * r);
                    uTot -= qiqj * (1 - ec) / r;
                }
                double s2 = 0, s4 = 0, a2 = 0, a4 = 0, alpha66 = 0, e = 0;
                if (alpha6 > 0) {
                    s2 = 1 / r2;
                    s4 = s2 * s2;
                    double alpha62 = alpha6 * alpha6;
                    a2 = r2 * alpha62;
                    a4 = a2 * a2;
                    alpha66 = alpha62 * alpha62 * alpha62;
                    e = Math.exp(-a2);
                    double ureal = -Bij * alpha66 * (1 + a2 + a4 / 2) * e / (a4 * a2);
                    double ufull = -Bij * s4 * s2;
                    uTot += -ufull + ureal;
                }
                if (!bondingInfo.isOnlyRigidMolecules()) {
                    // skip any non-bonded pair
                    if (!bondingInfo.skipBondedPair(false, atom, atoms.get(jAtom))) continue;
                    double du = 0;
                    if (alpha > 0)
                        du = -qiqj * (twoosqrtpi * Math.exp(-alpha * alpha * r2) * alpha + (ec - 1) / r) / r2;
                    if (alpha6 > 0) {
                        double a6 = a4 * a2;
                        double dureal = Bij * alpha66 * (6 + 6 * a2 + 3 * a4 + a6) * e / a6;
                        double dufull = 6 * Bij * s2 * s4;
                        du += -dufull + dureal;
                    }
                    if (doForces) {
                        forces[iAtom].PEa1Tv1(du, dr);
                        forces[jAtom].PEa1Tv1(-du, dr);
                    }
                    virialTot += du * r2;
                }
            }
        }
        return uTot;
    }

    @Override
    public double computeOneOld(IAtom atom) {
        return computeOneInternal(atom, true);
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        return 0;
    }

    @Override
    public double computeOne(IAtom atom) {
        return computeOneInternal(atom, false);
    }

    private double computeOneInternal(IAtom atom, boolean oldEnergy) {
        int iAtom = atom.getLeafIndex();
        int iType = atom.getType().getIndex();
        double qi = chargesByType[iType];
        if (qi==0 && B6[iType][iType]==0) return 0;

        double u = 0;
        double q2Sum = qi*qi;
        u -= alpha/SQRT_PI*q2Sum;
        double vol = box.getBoundary().volume();
        if (alpha6 > 0 && B6[iType][iType]!=0) {
            double sumBii = B6[iType][iType];
            double sumBij = 0;
            for (int jType=0; jType<atomCountByType.length; jType++) {
                sumBij += atomCountByType[jType]*B6[iType][jType];
            }
            double alpha63 = alpha6 * alpha6 * alpha6;
            u -= SQRT_PI * PI/(6*vol) * alpha63 * sumBij;
            u += (1.0/12.0) * alpha63 * alpha63 * sumBii;
        }
        Vector bs = box.getBoundary().getBoxSize();int kxMax = (int)(0.5 * bs.getX(0) / PI * kCut);
        // cube instead of sphere, so conservatively big
        int[] kMax = {
                kxMax,
                (int)(0.5 * bs.getX(1) / PI * kCut),
                (int)(0.5 * bs.getX(2) / PI * kCut)
        };
        int[] nk = {kMax[0] + 1, 2*kMax[1] + 1, 2*kMax[2] + 1};
        // we could save eik too... would need to update after accepted move, revert after rejection
        for (int a=0; a<3; a++) {
            double fac = 2.0*PI/bs.getX(a);
            if (a==0) q2Sum += qi*qi;
            int idx = iAtom*nk[a];
            if (a>0) idx += kMax[a];
            Vector ri = atom.getPosition();
            eik[a][idx] = 1;
            Complex.ONE.intoArray(eik[a], idx);
            new Complex(Math.cos(fac * ri.getX(a)), Math.sin(fac * ri.getX(a)))
                    .intoArray(eik[a], idx + 1);
            for (int i=2; i<=kMax[a]; i++) {
                eik[a][idx+i] = eik[a][idx+1] * eik[a][idx+i-1];
            }
            if (a==0) continue;
            for (int i=1; i<=kMax[a]; i++) {
                Complex.fromArray(eik[a], idx+i).conjugate().intoArray(eik[a], idx - 1);
            }
        }

        double fourierSum = 0;
        double fourierSum6 = 0;
        for (int ik = 0; ik < this.ik.size(); ik++) {
            int ikx = this.ik.getInt(ik * 3);
            int iky = this.ik.getInt(ik * 3 + 1);
            int ikz = this.ik.getInt(ik * 3 + 2);

            if (!oldEnergy) {
                Complex sFacMinus = Complex.fromArray(sFac, ik).minus(Complex.fromArray(dsFac, ik));
                // energy without this atom
                fourierSum -= fExp[ik] * (sFacMinus.times(sFacMinus.conjugate())).real();
                Complex.fromArray(dsFac, ik).changeSign().intoArray(dsFac, ik);
                if (alpha6 > 0) {
                    for (int kB=0; kB<=6; kB++) {
                        sFacMinus = Complex.fromArray(sFacB[kB], ik).minus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacMinus2 = Complex.fromArray(sFacB[6-kB], ik).minus(Complex.fromArray(dsFacB[6-kB], ik));
                        // energy without this atom
                        fourierSum6 -= f6Exp[ik] * (sFacMinus.times(sFacMinus2.conjugate())).real();
                    }
                    for (int kB=0; kB<=6; kB++) {
                        Complex.fromArray(dsFacB[kB], ik).changeSign().intoArray(dsFacB[kB], ik);
                    }
                }
            }
            Complex iContrib = Complex.fromArray(eik[0], iAtom * nk[0] + ikx)
                    .times(Complex.fromArray(eik[1], iAtom * nk[1] + kMax[1] + iky))
                    .times(Complex.fromArray(eik[2], iAtom * nk[2] + kMax[2] + ikz));
            Complex.fromArray(dsFac, ik).plus(iContrib.times(qi)).intoArray(dsFac, ik);
            if (alpha6 > 0) {
                for (int kB = 0; kB <= 6; kB++) {
                    Complex.fromArray(dsFacB[kB], ik).plus(iContrib.times(b6[iType][kB]))
                            .intoArray(dsFacB[kB], ik);
                }
            }

            if (oldEnergy) {
                Complex sFacMinus = Complex.fromArray(sFac, ik).minus(Complex.fromArray(dsFac, ik));
                fourierSum += fExp[ik] *
                        (Complex.fromArray(sFac, ik).times(Complex.fromArray(sFac, ik).conjugate()).real()
                        - sFacMinus.times(sFacMinus.conjugate()).real());

                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 6; kB++) {
                        sFacMinus = Complex.fromArray(sFacB[kB], ik).minus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacMinus2 = Complex.fromArray(sFacB[6-kB], ik).minus(Complex.fromArray(dsFacB[6-kB], ik));
                        fourierSum6 += f6Exp[ik] *
                                (Complex.fromArray(sFacB[kB], ik).times(Complex.fromArray(sFacB[kB], ik).conjugate()).real()
                                        - sFacMinus.times(sFacMinus2.conjugate()).real());
                    }
                }
            } else {
                Complex sFacNew = Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik));
                fourierSum += fExp[ik] * sFacNew.times(sFacNew.conjugate()).real();
                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 6; kB++) {
                        sFacNew = Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik));
                        Complex sFacNew2 = Complex.fromArray(sFacB[6-kB], ik).plus(Complex.fromArray(dsFacB[6-kB], ik));
                        fourierSum6 += f6Exp[ik] * sFacNew.times(sFacNew2.conjugate()).real();
                    }
                }
            }
        }
        u += fourierSum;
        u += fourierSum6;
        return u;
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        return 0;
    }

    @Override
    public void processAtomU(double fac) {
        for (int ik = 0; ik < nWaveVectors; ik++) {
            if (fac == 1) {
                Complex.fromArray(sFac, ik).plus(Complex.fromArray(dsFac, ik)).intoArray(sFac, ik);
                if (alpha6 > 0) {
                    for (int kB = 0; kB <= 6; kB++) {
                        Complex.fromArray(sFacB[kB], ik).plus(Complex.fromArray(dsFacB[kB], ik)).intoArray(sFac, ik);
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
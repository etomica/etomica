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
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import static etomica.math.SpecialFunctions.erfc;
import static etomica.math.SpecialFunctions.factorial;
import static java.lang.Math.PI;

public class PotentialComputeEwaldFourier implements PotentialCompute {
    private static final double SQRT_PI = Math.sqrt(PI);
    private final Box box;
    protected double[] uAtom;
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
    private Complex[] sFacAtom = new Complex[0]; // Complex for each atom
    private Complex[] sFac = new Complex[0]; // Complex for each kVector
    private final Complex[][] sFacB = new Complex[7][0]; // 7 arrays of, Complex for each kVector

    // Array for each spacial dimension, then flattened array of (num kVectors in that dimension)*Complex for each atom
    private final Complex[][] eik = new Complex[3][0];

    private Complex[] dsFac = new Complex[0]; // Complex for each kVector
    private final Complex[][] dsFacB = new Complex[7][0]; // 7 arrays of, Complex for each kVector
    private double[] fExp; // double for each kVector
    private double[] f6Exp; // double for each kVector
    private int nWaveVectors;


    private void setArraySizes(int numAtoms, int numKVectors, int[] dimKVectors) {
        if (numAtoms > sFacAtom.length) {
            sFacAtom = new Complex[numAtoms];
            Arrays.setAll(sFacAtom, i -> new Complex());

            forces = new Vector[numAtoms];
            Arrays.setAll(forces, i -> space.makeVector());
        }

        for (int i = 0; i < eik.length; i++) {
            if (dimKVectors[i] * numAtoms > eik[i].length) {
                eik[i] = new Complex[numAtoms * dimKVectors[i]];
                Arrays.setAll(eik[i], j -> new Complex());
            }
        }

        if (numKVectors > sFac.length) {
            sFac = new Complex[numKVectors];
            Arrays.setAll(sFac, i -> new Complex());

            for (int j = 0; j < sFacB.length; j++) {
                sFacB[j] = new Complex[numKVectors];
                Arrays.setAll(sFacB[j], i -> new Complex());
            }

            dsFac = new Complex[numKVectors];
            Arrays.setAll(dsFac, i -> new Complex());

            for (int j = 0; j < dsFacB.length; j++) {
                dsFacB[j] = new Complex[numKVectors];
                Arrays.setAll(dsFacB[j], i -> new Complex());
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
        return 0;
    }

    @Override
    public void updateAtom(IAtom atom) {

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

        // TODO intramolecular

        double kCut2 = kCut * kCut;
        Vector bs = box.getBoundary().getBoxSize();
        int kxMax = (int)(0.5 * bs.getX(0) / PI * kCut);
        // cube instead of sphere, so conservatively big
        int[] kMax = {
                kxMax,
                (int)(0.5 * bs.getX(1) / PI * kCut),
                (int)(0.5 * bs.getX(2) / PI * kCut)
        };
        int[] nk = {kMax[0] + 1, 2*kMax[1] + 1, 2*kMax[2] + 1};
        int nkTot = ((2*kxMax + 1) * nk[1] * nk[2] - 1) / 2;

        this.setArraySizes(numAtoms, nkTot, nk);
        Arrays.stream(forces).forEach(v -> v.E(0));

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
                eik[a][idx] = new Complex(1, 0);
                if (nk[a]==1) continue;
                int iType = atoms.get(iAtom).getType().getIndex();
                if (chargesByType[iType] == 0 && B6[iType][iType] == 0) continue;
                Vector ri = atoms.get(iAtom).getPosition();
                eik[a][idx+1] = new Complex(Math.cos(fac*ri.getX(a)), Math.sin(fac*ri.getX(a)));
                for (int i=2; i<=kMax[a]; i++) {
                    eik[a][idx+i] = eik[a][idx+1].times(eik[a][idx+i-1]);
                }
                if (a==0) continue;
                for (int i=1; i<=kMax[a]; i++) {
                    eik[a][idx-i] = eik[a][idx+i].conjugate();
                }
            }
        }

        double coeff = 4 * PI / vol;
        double coeffB2 = -2.0* SQRT_PI * PI * alpha6*alpha6*alpha6/(3.0*vol);
        double fourierSum = 0;
        double fourierSum6 = 0;
        double virialSum = 0;
        double virialSum6 = 0;
        int ik = 0;
        kBasis.E(2 * PI);
        kBasis.DE(bs);

        for (int ikx=0; ikx<=kxMax; ikx++) {
            double kx = ikx*kBasis.getX(0);
            double kx2 = kx*kx;
            double kyCut2 = kCut2 - kx2;
            boolean xpositive = ikx>0;
            int kyMax = (int)(0.5*bs.getX(1)*Math.sqrt(kyCut2)/PI);
            for (int iky=-kyMax; iky<=kyMax; iky++) {
                if (!xpositive && iky<0) continue;
                boolean ypositive = iky>0;
                double ky = iky*kBasis.getX(1);
                double kxy2 = kx2 + ky*ky;
                int kzMax = (int)(0.5*bs.getX(2)*Math.sqrt(kCut2 - kxy2)/PI);
                for (int ikz=-kzMax; ikz<=kzMax; ikz++) {
                    if (!xpositive && !ypositive && ikz<=0) continue;
                    double kz = ikz*kBasis.getX(2);
                    double kxyz2 = kxy2 + kz*kz;
                    double expthing = Math.exp(-0.25*kxyz2/(alpha*alpha));
                    sFac[ik] = Complex.ZERO;
                    for (int kB=0; kB<=6; kB++) sFacB[kB][ik] = Complex.ZERO;
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
                    for (int iAtom=0; iAtom<numAtoms; iAtom++) {
                        int iType = atoms.get(iAtom).getType().getIndex();
                        double qi = chargesByType[iType];
                        double Bii = B6[iType][iType];
                        if (qi==0 && Bii == 0) {
                            sFacAtom[iAtom] = Complex.ZERO;
                            continue;
                        }
                        Complex iContrib = eik[0][iAtom*nk[0]+ikx]
                                .times(eik[1][iAtom*nk[1]+kMax[1]+iky])
                                .times(eik[2][iAtom*nk[2]+kMax[2]+ikz]);
                        sFacAtom[iAtom] = iContrib;
                        sFac[ik] = sFac[ik].plus(iContrib.times(qi));
                        if (alpha6 > 0) {
                            for (int kB=0; kB<=6; kB++) {
                                sFacB[kB][ik] = sFacB[kB][ik].plus(iContrib.times(b6[iType][kB]));
                            }
                        }
                    }
                    double x = (sFac[ik].times(sFac[ik].conjugate())).real();
                    fourierSum += fExp[ik] * x;
                    double kdfdk = 0;
                    double dfqdV = 0, df6dV = 0;
                    if (alpha>0) {
                        kdfdk = -(2 + kxyz2/(2*alpha*alpha))*fExp[ik];
                        dfqdV = -fExp[ik]/vol - kdfdk/(3*vol);
                        virialSum += 3*vol*dfqdV*x;
                    }
                    if (alpha6>0) {
                        df6dV = -f6Exp[ik]/vol - hdf6dh/(3*vol);
                        for (int kB=0; kB<=6; kB++) {
                            double y = sFacB[kB][ik].times(sFacB[6 - kB][ik].conjugate()).real();
                            fourierSum6 += f6Exp[ik] * y;
                            virialSum6 += 3*vol*df6dV*y;
                        }
                    }
                    if (doForces) {
                        for (int iAtom=0; iAtom<numAtoms; iAtom++) {
                            int iType = atoms.get(iAtom).getType().getIndex();
                            double coeffki = alpha==0 ? 0 :
                                    2*fExp[ik]*chargesByType[iType]*(sFacAtom[iAtom].times(sFac[ik].conjugate()).imaginary());
                            double coeffki6 = 0;
                            if (alpha6 > 0) {
                                for (int kB=0; kB<=6; kB++) {
                                    coeffki6 += 2*f6Exp[ik]*b6[iType][kB]*(sFacAtom[iAtom].times(sFacB[6-kB][ik].conjugate())).imaginary();
                                }
                            }
                            coeffki += coeffki6;
                            forces[iAtom].PE(Vector.of(coeffki * kx, coeffki * ky, coeffki * kz));
                        }
                    }
                    ik++;
                }
            }
        }

        this.nWaveVectors = ik;
        uTot += fourierSum + fourierSum6;
        virialTot += virialSum - virialSum6;

        return uTot;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return 0;
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        return 0;
    }

    @Override
    public double computeOne(IAtom iAtom) {
        return 0;
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        return 0;
    }

    @Override
    public void processAtomU(double fac) {

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

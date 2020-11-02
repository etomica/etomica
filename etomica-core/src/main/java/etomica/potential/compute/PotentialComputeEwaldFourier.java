package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
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

import static etomica.math.SpecialFunctions.factorial;

public class PotentialComputeEwaldFourier implements PotentialCompute {
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
    private Complex[] sFacAtom; // Complex for each atom
    private Complex[] sFac; // Complex for each kVector
    private final Complex[][] sFacB = new Complex[7][]; // 7 arrays of, Complex for each kVector

    // Array for each spacial dimension, then flattened array of (num kVectors in that dimension)*Complex for each atom
    private final Complex[][] eik = new Complex[3][];

    private Complex[] dsFac; // Complex for each kVector
    private final Complex[][] dsFacB = new Complex[7][]; // 7 arrays of, Complex for each kVector
    private double[] fExp; // double for each kVector
    private double[] f6Exp; // double for each kVector


    private void setArraySizes(int numAtoms, int numKVectors, int[] dimKVectors) {
        if (numAtoms > sFacAtom.length) {
            sFacAtom = new Complex[numAtoms];
            Arrays.setAll(sFacAtom, i -> new Complex());
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
        return new Vector[0];
    }

    @Override
    public double getLastVirial() {
        return 0;
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
        return 0;
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
        return null;
    }
}

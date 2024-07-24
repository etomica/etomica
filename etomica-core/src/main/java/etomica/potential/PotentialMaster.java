/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;

public class PotentialMaster implements PotentialCompute {
    protected final BondingInfo bondingInfo;
    protected final IPotential2[][] pairPotentials;
    protected double[][][] pScale;
    protected final Box box;
    protected final Vector dr;
    protected double u;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected Vector zero;
    protected double virialTot = Double.NaN, energyTot = Double.NaN;
    protected Vector[] forces;
    protected final Space space;
    protected double minR2;

    protected final boolean isPureAtoms;
    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;
    protected final double[] reduction;
    boolean first = true;
    public long tAll, numAll;
    public long numMC, tMC;


    public PotentialMaster(SpeciesManager sm, Box box, BondingInfo bondingInfo) {
        this(sm, box, bondingInfo, new IPotential2[sm.getAtomTypeCount()][sm.getAtomTypeCount()]);
    }
    public PotentialMaster(SpeciesManager sm, Box box, BondingInfo bondingInfo, IPotential2[][] pairPotentials) {
        space = box.getSpace();
        this.bondingInfo = bondingInfo;
        this.pairPotentials = pairPotentials;
        pScale = new double[sm.getAtomTypeCount()][sm.getAtomTypeCount()][];
        reduction = new double[pairPotentials.length];
        this.box = box;
        dr = box.getSpace().makeVector();

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        forces = new Vector[0];

        isPureAtoms = sm.isPureAtoms();

        this.atomCountByType = new int[sm.getAtomTypeCount()];
        for (ISpecies s : sm.getSpeciesArray()) {
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

                int newAtoms = e.getMolecule().getType().getLeafAtomCount();
                int nowAtoms = box.getLeafList().size();
                if (nowAtoms > uAtom.length) {
                    double[] uAtomNew = new double[nowAtoms];
                    System.arraycopy(uAtom, 0, uAtomNew, 0, nowAtoms - newAtoms);
                    uAtom = uAtomNew;
                } else {
                    Arrays.fill(uAtom, nowAtoms - newAtoms, nowAtoms, 0);
                }
            }

            @Override
            public void boxNumberMolecules(BoxMoleculeCountEvent e) {
                int n = e.getCount();

                int nowAtoms = box.getLeafList().size();
                int newAtoms = e.getSpecies().getLeafAtomCount() * n;
                if (nowAtoms + newAtoms > uAtom.length) {
                    double[] uAtomNew = new double[nowAtoms + newAtoms];
                    System.arraycopy(uAtom, 0, uAtomNew, 0, nowAtoms);
                    uAtom = uAtomNew;
                }
            }

            @Override
            public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) {
                int oldIndex = e.getIndex();
                int newIndex = e.getAtom().getLeafIndex();
                uAtom[newIndex] = uAtom[oldIndex];
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]--;
                }
            }
        });
    }

    @Override
    public boolean needForcesForVirial() {
        return !isPureAtoms;
    }

    @Override
    public void init() {
        Vector L = box.getBoundary().getBoxSize();
        double minL = Double.POSITIVE_INFINITY;
        for (int i=0; i<L.getD(); i++) {
            minL = Math.min(minL, L.getX(0));
        }
        minR2 = 0.25*minL*minL;
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
        return energyTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, IPotential2 p12) {
        setPairPotential(atomType1, atomType2, p12, new double[]{1,0,0,0});
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, IPotential2 p12, double[] ijScale){
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;
        pScale[atomType2.getIndex()][atomType1.getIndex()]= ijScale;
        pScale[atomType1.getIndex()][atomType2.getIndex()]= ijScale;
    }

    public IPotential2[][] getPairPotentials() {
        return pairPotentials;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    protected double handleComputeAll(boolean doForces, IAtom iAtom, IAtom jAtom, int iAtomNum, int jAtomNum, Vector ri, Vector rj, Vector jbo, IPotential2 ip2, PotentialCallback pc){
        double u =0;

        if (pc != null && pc.skipPair(iAtomNum, jAtomNum)) return 0;
        numAll++;
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double[] u012 = new double[3];
        double r2 = dr.squared();
        double multiplier = r2 < minR2 ? getpScaleMultiplier(iAtom, jAtom) : 1;
        ip2.u012add(r2, u012);
        double uij = multiplier*u012[0];
        u += uij;
        uAtom[iAtomNum] += 0.5 * u;
        uAtom[jAtomNum] += 0.5 * u;
        double duij;
        duij = multiplier * u012[1];
        virialTot += duij;
        if (doForces && duij != 0) {
            dr.TE(duij / r2);
            forces[iAtomNum].PE(dr);
            forces[jAtomNum].ME(dr);
            if(forces[iAtomNum].equals(Double.NaN) || forces[jAtomNum].equals(Double.NaN)){

                throw new RuntimeException("oops ");

            }
        }
       // System.out.println(iAtom + " " + iAtom.getPosition() + jAtom + " " + jAtom.getPosition() + " "+ u);
        return u;
    }


    protected final void zeroArrays(boolean doForces) {
        virialTot = 0;

        int numAtoms = box.getLeafList().size();
        if (doForces && numAtoms > forces.length) {
            int oldLength = forces.length;
            forces = Arrays.copyOf(forces, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) forces[i] = box.getSpace().makeVector();
        }
        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
        }
        for (int i = 0; i < numAtoms; i++) {
            uAtom[i] = 0;
            if (doForces) forces[i].E(0);
        }
    }


    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        double[] uAtomOld = new double[uAtom.length];
        boolean debug = false;
        if (debug) System.arraycopy(uAtom, 0, uAtomOld, 0, uAtom.length);
        zeroArrays(doForces);
        IAtomList atoms = box.getLeafList();
        double u = 0;
        Boundary boundary = box.getBoundary();
        long t1 = System.nanoTime();

        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            int iType = iAtom.getType().getIndex();
            IPotential2[] ip1 = pairPotentials[iType];
            for (int j = i + 1; j < atoms.size(); j++) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Vector rj = jAtom.getPosition();
                IPotential2 ip2 = ip1[jType];
                if (ip2 == null) continue;

                if (bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom)) continue;

                dr.Ev1Mv2(jAtom.getPosition(), ri);
                boundary.nearestImage(dr);
                u += handleComputeAll(doForces, iAtom, jAtom, i, j, zero, dr, zero, ip2, pc);
            }
        }
        tAll += System.nanoTime() - t1;

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        u += uCorrection[0];
        virialTot += duCorrection[0];

        if (doForces && !isPureAtoms) {
            virialTot += PotentialCompute.computeVirialIntramolecular(forces, box);
        }

        if (debug && uAtom.length == uAtomOld.length && !first) {
            boolean success = true;
            for (int i = 0; i < uAtom.length; i++) {
                if (Math.abs(uAtom[i] - uAtomOld[i]) > 1e-9) {
                    success = false;
                }
            }
            if (!success) throw new RuntimeException("oops");
        }
        first = false;
        energyTot = u;
        return u;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u + uAtom[iAtom.getLeafIndex()] * 2;
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        double u = 0;
        for (IAtom atom : atoms) {
            u += uAtom[atom.getLeafIndex()] * 2;
            u += this.computeOneTruncationCorrection(atom.getLeafIndex());
        }
        u -= computeIntraAtoms(atoms);
        return u;
    }

    protected boolean arrayContains(IAtom a, int startIdx, IAtom... atoms) {
        for (int i=startIdx; i<atoms.length; i++) {
            if (atoms[i] == a) return true;
        }
        return false;
    }

    protected double computeIntraAtoms(IAtom... atoms) {
        // XXX when using neighbor lists, we should only include neighbors
        //     or, computeOneInternal needs to avoid double-counting
        double uIntra = 0;
        for (int i = 0; i < atoms.length; i++) {
            IAtom atom1 = atoms[i];
            for (int j = i + 1; j < atoms.length; j++) {
                IAtom atom2 = atoms[j];
                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) continue;
                IPotential2 pij = pairPotentials[atom1.getType().getIndex()][atom2.getType().getIndex()];
                if (pij == null) continue;
                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                box.getBoundary().nearestImage(rij);
                double uij = pij.u(rij.squared());
                uIntra += uij;
            }
        }
        return uIntra;
    }
    protected double getpScaleMultiplier(IAtom iAtom, IAtom jAtom){
        double[][] pScale1 = pScale[iAtom.getType().getIndex()];
        int n = bondingInfo.n(false, iAtom, jAtom);
        if(n>3){
            n = 3;
        }
        int nStar = jAtom.getType().getIndex();
        return pScale1[nStar][n];
    }
    protected double handleComputeOne( IAtom iAtom, IAtom jAtom, Vector ri, Vector rj, Vector jbo, IPotential2 ip2){
        numMC++;
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double r2 = dr.squared();
        //System.out.println("Non");
        double uij = ip2.u(r2);
        if (uij == 0) {
            return 0;
        }
        double multiplier = r2 < minR2 ? getpScaleMultiplier(iAtom, jAtom) : 1;
        uij *= multiplier;
        uAtomsChanged.add(jAtom.getLeafIndex());

        if (duAtomMulti) {
            duAtom.plusEquals(iAtom.getLeafIndex(), 0.5 * uij);
            duAtom.plusEquals(jAtom.getLeafIndex(), 0.5 * uij);
        } else {
            duAtom.plusEquals(0, 0.5 * uij);
            duAtom.add(0.5 * uij);
        }
       //  System.out.println("Non " + uij);
        //System.out.println(iAtom.getType() + " "+ iAtom.getPosition()+" " + jAtom.getType() + " "+jAtom.getPosition()+" " + uij +" "+  Math.exp(-uij/291));
        return uij;

    }

    @Override
    public double computeOne(IAtom iAtom) {
        this.duAtomMulti = false;
        int i = iAtom.getLeafIndex();
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
        uAtomsChanged.ensureCapacity(atoms.size());
        uAtomsChanged.add(i);
        duAtom.clear();
        duAtom.ensureCapacity(atoms.size());
        duAtom.add(0);

        double u = computeOneInternal(iAtom);
        u += this.computeOneTruncationCorrection(i);
        return u;
    }

    protected double computeOneInternal(IAtom atom) {
        return computeOneInternal(atom, 0, new IAtom[0]);
    }

    protected double computeOneInternal(IAtom atom, int startExcludeIdx, IAtom... excludedAtoms) {
        int iType = atom.getType().getIndex();
        int i = atom.getLeafIndex();
        IPotential2[] ip = pairPotentials[iType];
        double u = 0;
        Boundary boundary = box.getBoundary();
        long t1 = System.nanoTime();
        for (int j = 0; j < box.getLeafList().size(); j++) {
            if (i == j) continue;
            IAtom jAtom = box.getLeafList().get(j);
            int jType = jAtom.getType().getIndex();
            IPotential2 pij = ip[jType];
            if (pij == null) continue;
            if (arrayContains(jAtom, startExcludeIdx, excludedAtoms)) continue;
            if (bondingInfo.skipBondedPair(isPureAtoms, atom, jAtom)) continue;
            dr.Ev1Mv2(jAtom.getPosition(), atom.getPosition());
            boundary.nearestImage(dr);
            u += handleComputeOne(atom, jAtom, zero, dr,  zero,pij);
        }
        tMC += System.nanoTime() - t1;
        return u;
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        duAtomMulti = true;
        uAtomsChanged.clear();
        duAtom.setAll(0);
        duAtom.ensureCapacity(box.getLeafList().size());
        double u = 0;

        for (int i=0; i<atoms.length; i++) {
            IAtom atom = atoms[i];
            uAtomsChanged.add(atom.getLeafIndex());
            u += computeOneInternal(atom, i+1, atoms);
            u += computeOneTruncationCorrection(atom.getLeafIndex());
        }
        if (u == Double.POSITIVE_INFINITY) return u;
       // System.out.println("Uexcess "+ u);
        return u;
    }

    @Override
    public void processAtomU(double fac) {
        for (int j = 0; j < uAtomsChanged.size(); j++) {
            int jj = uAtomsChanged.getInt(j);
            if (duAtomMulti) {
                uAtom[jj] += fac * duAtom.getDouble(jj);
                duAtom.replace(jj, 0);
            } else {
                uAtom[jj] += fac * duAtom.getDouble(j);
            }
        }
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {
                init();
            }

            @Override
            public void integratorStepStarted(IntegratorEvent e) {

            }

            @Override
            public void integratorStepFinished(IntegratorEvent e) {

            }
        };
    }

    public void computeAllTruncationCorrection(double[] uCorrection, double[] duCorrection) {
        if (!doAllTruncationCorrection) {
            return;
        }
        double uCor = 0;
        double duCor = 0;
        for (int i = 0; i < atomCountByType.length; i++) {
            for (int j = i; j < atomCountByType.length; j++) {
                IPotential2 p = pairPotentials[i][j];
                if (p == null) continue;
                int numPairs;
                if (j == i) {
                    numPairs = atomCountByType[i] * (atomCountByType[j] - 1) / 2;
                } else {
                    numPairs = atomCountByType[i] * atomCountByType[j];
                }
                double pairDensity = numPairs / box.getBoundary().volume();

                double[] u = new double[1];
                double[] du = new double[1];
                p.u01TruncationCorrection(space, u, du);
                uCor += pairDensity * u[0];
                duCor += pairDensity * du[0];

            }
        }
        uCorrection[0] = uCor;
        duCorrection[0] = duCor;
    }

    public double computeOneTruncationCorrection(int iAtom) {
        double integral =0;
        if (!doOneTruncationCorrection) {
            return 0;
        }
        int iType = box.getLeafList().get(iAtom).getType().getIndex();
        double uCorrection = 0;
        for (int j = 0; j < atomCountByType.length; j++) {
            IPotential2 p = pairPotentials[iType][j];
            double pairDensity;
            if (iType == j) {
                pairDensity = (atomCountByType[j] - 1) / box.getBoundary().volume();
            } else {
                pairDensity = atomCountByType[j] / box.getBoundary().volume();
            }
            if(p == null){
                integral = 0;
            }else {
                integral = p.integral(space, p.getRange());
            }

            uCorrection += pairDensity * integral;
        }
        return uCorrection;
    }

}

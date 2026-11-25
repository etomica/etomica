/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialMolecular;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.virial.BoxCluster;

/**
 * Molecular pair potential for Mayer sampling of large flexible molecules (star polymers), an
 * alternative to {@link etomica.virial.PotentialMoleculePairCell} that updates energies
 * incrementally instead of recomputing every pair after each move.
 * <p>
 * Like PotentialMoleculePairCell, it computes all the energies for the box at once (the
 * intramolecular energy of each molecule and the interaction energy of each pair of molecules)
 * whenever the box's coordinate-pair ID changes, and returns stored values while the ID stays the
 * same.  It keeps the last numStates calculations (positions and energies); when the ID goes back
 * to one of them (a rejected trial), that calculation is reused.  A new calculation starts from the
 * stored ones, separately for each molecule:
 * <ul>
 * <li>It compares the positions with those of each stored calculation to find which atoms moved.
 * Each stored calculation has correct energies for its own positions, so any of them can be the
 * starting point.  This matters because the potential is not always called after a rejection
 * (ClusterSum keeps its own cache), so the most recent calculation may be a rejected trial.</li>
 * <li>If no atom of the molecule moved relative to a stored calculation, its intramolecular energy
 * is reused.</li>
 * <li>If more than minRigidTest of its atoms moved (for example all of them), it is reused when
 * the molecule moved rigidly, which is tested exactly: all distances from every atom to four non-coplanar reference atoms must be unchanged
 * (to round-off).  Positions are fixed by their distances to four non-coplanar points, so this
 * holds only for a rigid translation and/or rotation.  This also covers a molecule whose atoms
 * differ only by round-off, as after a rejected rotation that the move undoes by rotating back.</li>
 * <li>Otherwise, for the stored calculation with the fewest moved atoms (set M), the intramolecular
 * energy changes by the sum, over pairs with at least one atom in M, of the new minus the old pair
 * energy.  This is exact however M moved.  Pairs of a moved atom with unmoved non-core atoms are
 * found with a grid of the unmoved atoms (which are at the same positions before and after), so
 * the cost scales with the number of moved atoms.  If every atom moved non-rigidly relative to every stored
 * calculation, the molecule's intramolecular energy is recomputed in full.  It is also recomputed
 * in full when the update would lose precision: when the stored energy or the summed pair terms
 * are more than 1e4 times the result (large terms cancelling, from overlapping atoms in the
 * stored or new configuration), or when anything in the update is infinite or NaN (coincident
 * atoms, as in the starting configuration of the Mayer-sampling drivers).</li>
 * <li>The interaction energy of two molecules is reused when neither moved relative to a stored
 * calculation.  When one is unchanged and only some atoms of the other moved (an angle move), it
 * changes by the sum over the moved atoms of their new minus old energy with the unchanged
 * molecule.  Otherwise it is recomputed in full: one molecule's non-core atoms are placed in a grid
 * (cell size at least the cutoff) covering only its bounding box, so there is no fixed box, and
 * atoms of the other molecule outside that box (plus the cutoff) are skipped at once.  It is also
 * recomputed in full when an update would lose precision (the interaction energy is often huge,
 * from overlapping cores).</li>
 * <li>Atoms count as moved only if a coordinate changed by more than MOVE_TOLERANCE * (1 + |x|), so
 * that positions restored to round-off after a rejected rotation count as unchanged.</li>
 * </ul>
 * Atoms of the core type interact with every atom (their potentials may have infinite range).
 * All other pairs must have a finite range (getRange), which is used as the cutoff.  Distances are
 * plain differences (no periodic images), as needed for isolated molecules.
 * <p>
 * <b>Stored calculations and ClusterSum caching.</b>  Updates are cheapest when they start from
 * the accepted configuration.  ClusterSum caches its value by coordinate-pair ID, and with that
 * caching on it restores its value after a rejected trial without calling the potential, so the
 * potential does not learn that the trial was rejected: after k rejections in a row the accepted
 * configuration is k+1 calculations back.  Two ways to handle this:
 * <ul>
 * <li>Turn ClusterSum's caching off ({@code setCaching(false)}).  The moves ask for the cluster
 * value (or the intramolecular energy) at the start of each trial, before moving anything, so the
 * potential is then called at the restored ID after every rejection and finds the accepted
 * configuration among its stored ones.  Two stored calculations (setNumStates(2)) are then
 * enough.  Recomputing the cluster value from the stored energies is cheap.  This was the fastest
 * setup for 8-arm, 60-monomer stars (B2): about 7% faster than four stored with caching on.</li>
 * <li>Keep the caching and store more calculations (the default, 4), so the accepted
 * configuration usually survives a few rejections in a row.  With caching on and only two stored,
 * updates often start from a rejected trial, which was about 25% slower.</li>
 * </ul>
 * Either way the energies are correct: any stored calculation is a valid starting point, and the
 * choice only affects how much work an update takes.  A move that does not ask for the energy
 * before moving (for example with skipW set) only makes updates more expensive.
 * <p>
 * To limit round-off drift, a full recalculation is done every fullInterval calculations.  With
 * debugCount > 0, each of the next debugCount incremental calculations is compared with a full
 * recalculation, and a mismatch throws an exception.
 */
public class PotentialMoleculePairStar implements IPotentialMolecular {

    protected final IPotential2[][] atomPotentials;
    protected final AtomType coreType;
    protected BoxCluster box;
    protected int nAtoms, nMol;
    protected int[][] molCores, molOthers;   // leaf indices of core and non-core atoms of each molecule
    protected int[] atomMol, atomType;
    protected boolean[] isCore;
    protected double rc, rc2;
    protected State[] states;
    protected State cur;
    protected int numStates = 4;
    protected boolean[] moved;
    protected int[] movedList;
    // grid for the interaction energies
    protected int[] cellHead = new int[0], cellNext, unmovedList;
    protected final double[] gMin = new double[3], gMax = new double[3], gCellInv = new double[3];
    protected final int[] gNc = new int[3];
    protected boolean gridEmpty;
    protected double gridAbs;
    protected final int maxCellsPerDim = 64;

    public int fullInterval = 1000;
    /** rigidity is tested relative to a stored calculation when more than this many atoms moved */
    public int minRigidTest = 8;
    protected boolean forceFull;
    public long debugCount = 0;
    public double debugTolerance = 1e-8;
    /** counts of full and incremental calculations, and of restores after rejected trials */
    public long numFull, numIncremental, numRestored;
    /** counts of how each molecule's intramolecular energy was obtained in incremental calculations */
    public long numMolUnchanged, numMolRigid, numMolPartial, numMolFull;
    /** counts of how interaction energies were obtained in incremental calculations */
    public long numInterUnchanged, numInterPartial, numInterFull;
    /** total moved atoms over the partial updates (intramolecular and interaction) */
    public long sumMovedIntra, sumMovedInter;

    protected static class State {
        double[] r;
        double[] intra;
        double[][] inter;
        long id = -1;
        long lastUsed;

        State(int nAtoms, int nMol) {
            r = new double[3 * nAtoms];
            intra = new double[nMol];
            inter = new double[nMol][nMol];
        }
    }

    /**
     * Sets how many calculations are kept (at least 2; the default is 4).  Two is enough when
     * ClusterSum's caching is off; see the class comment.  Call before setBox.
     */
    public void setNumStates(int n) {
        numStates = Math.max(2, n);
    }

    public PotentialMoleculePairStar(SpeciesManager sm, AtomType coreType) {
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        atomPotentials = new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
        this.coreType = coreType;
    }

    public void setAtomPotential(AtomType atomType1, AtomType atomType2, IPotential2 p2) {
        atomPotentials[atomType1.getIndex()][atomType2.getIndex()] = p2;
        atomPotentials[atomType2.getIndex()][atomType1.getIndex()] = p2;
    }

    public IPotential2[][] getAtomPotentials() {
        return atomPotentials;
    }

    /**
     * Sets the box.  Call this after the molecules have been added and all atom potentials set.
     */
    public void setBox(BoxCluster box) {
        this.box = box;
        IMoleculeList molecules = box.getMoleculeList();
        IAtomList atoms = box.getLeafList();
        nMol = molecules.size();
        nAtoms = atoms.size();
        atomMol = new int[nAtoms];
        atomType = new int[nAtoms];
        isCore = new boolean[nAtoms];
        molCores = new int[nMol][];
        molOthers = new int[nMol][];
        rc = 0;
        for (int m = 0; m < nMol; m++) {
            IAtomList children = molecules.get(m).getChildList();
            int nc = 0;
            for (IAtom a : children) if (a.getType() == coreType) nc++;
            molCores[m] = new int[nc];
            molOthers[m] = new int[children.size() - nc];
            int ic = 0, io = 0;
            for (IAtom a : children) {
                int i = a.getLeafIndex();
                atomMol[i] = m;
                atomType[i] = a.getType().getIndex();
                isCore[i] = a.getType() == coreType;
                if (isCore[i]) molCores[m][ic++] = i;
                else molOthers[m][io++] = i;
            }
        }
        for (int i = 0; i < nAtoms; i++) {
            for (int j = 0; j < nAtoms; j++) {
                if (isCore[i] || isCore[j]) continue;
                IPotential2 p = atomPotentials[atomType[i]][atomType[j]];
                if (p == null) continue;
                if (Double.isInfinite(p.getRange())) {
                    throw new RuntimeException("non-core potentials need a finite range");
                }
                rc = Math.max(rc, p.getRange());
            }
        }
        rc2 = rc * rc;
        states = new State[numStates];
        for (int k = 0; k < numStates; k++) states[k] = new State(nAtoms, nMol);
        cur = null;
        moved = new boolean[nAtoms];
        movedList = new int[nAtoms];
        cellNext = new int[nAtoms];
        unmovedList = new int[nAtoms];
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {
        long thisCPairID = box.getCPairSet().getID();
        if (cur == null || thisCPairID != cur.id) {
            State found = null;
            for (State st : states) {
                if (st.id == thisCPairID) found = st;
            }
            if (found != null) {
                // back to a stored configuration (a rejected trial)
                cur = found;
                numRestored++;
            } else {
                // overwrite the least recently used calculation
                State next = states[0];
                for (State st : states) {
                    if (st.lastUsed < next.lastUsed) next = st;
                }
                next.id = -1;
                copyPositions(next.r);
                if (cur == null || forceFull || (numFull + numIncremental) % fullInterval == 0) {
                    computeFull(next);
                    forceFull = false;
                    numFull++;
                } else {
                    computeIncremental(next);
                    numIncremental++;
                    if (debugCount > 0) {
                        debugCount--;
                        State check = new State(nAtoms, nMol);
                        System.arraycopy(next.r, 0, check.r, 0, next.r.length);
                        computeFull(check);
                        compare(next, check);
                    }
                }
                next.id = thisCPairID;
                cur = next;
            }
            cur.lastUsed = numFull + numIncremental + numRestored;
        }
        int i1 = molecule1.getIndex(), i2 = molecule2.getIndex();
        return i1 == i2 ? cur.intra[i1] : cur.inter[i1][i2];
    }

    /**
     * Makes the next calculation (the next new coordinate-pair ID) a full recalculation.
     */
    public void forceFullNextCalculation() {
        forceFull = true;
    }

    /**
     * Recomputes every energy for the given state (whose positions are already set), from scratch.
     */
    public void computeFull(State s) {
        double[] r = s.r;
        for (int m = 0; m < nMol; m++) {
            double u = 0;
            int[] cores = molCores[m], others = molOthers[m];
            for (int a = 0; a < cores.length; a++) {
                int i = cores[a];
                for (int b = a + 1; b < cores.length; b++) u += uCore(r, i, cores[b]);
                for (int j : others) u += uCore(r, i, j);
            }
            for (int a = 0; a < others.length; a++) {
                int i = others[a];
                for (int b = a + 1; b < others.length; b++) u += uOther(r, i, others[b]);
            }
            s.intra[m] = u;
        }
        computeInter(s);
    }

    /** intramolecular energy of molecule m, from scratch */
    protected double intraFull(double[] r, int m) {
        double u = 0;
        int[] cores = molCores[m], others = molOthers[m];
        for (int a = 0; a < cores.length; a++) {
            int i = cores[a];
            for (int b = a + 1; b < cores.length; b++) u += uCore(r, i, cores[b]);
            for (int j : others) u += uCore(r, i, j);
        }
        for (int a = 0; a < others.length; a++) {
            int i = others[a];
            for (int b = a + 1; b < others.length; b++) u += uOther(r, i, others[b]);
        }
        return u;
    }

    /**
     * Computes the energies for state s (positions set) from the stored calculations; see the
     * class comment.
     */
    protected void computeIncremental(State s) {
        double[] rNew = s.r;
        for (int m = 0; m < nMol; m++) {
            int[] cores = molCores[m], others = molOthers[m];
            int molSize = cores.length + others.length;
            State unchanged = null, rigid = null, partial = null;
            int partialMoved = molSize;
            // most recently used first, so ties go to the most recent calculation
            State[] byRecency = states.clone();
            java.util.Arrays.sort(byRecency, (a, b) -> Long.compare(b.lastUsed, a.lastUsed));
            for (State c : byRecency) {
                if (c == s || c.id < 0) continue;
                int n = countMoved(c.r, rNew, cores, others);
                if (n == 0) {
                    unchanged = c;
                    break;
                }
                // Test for a rigid move whenever more than a few atoms moved, not only when all
                // did: a rejected rotation is undone by rotating back, which restores positions
                // only to round-off, leaving some (not all) atoms of an unmoved molecule changed.
                if (n > minRigidTest && rigid == null && isRigid(c.r, rNew, m)) {
                    rigid = c;
                    continue;
                }
                if (n == molSize) continue;
                if (n < partialMoved || (n == partialMoved && partial != null
                        && Math.abs(c.intra[m]) < Math.abs(partial.intra[m]))) {
                    // fewest moved atoms; among equals, the smallest energy (least cancellation)
                    partial = c;
                    partialMoved = n;
                }
            }
            if (unchanged != null) {
                s.intra[m] = unchanged.intra[m];
                numMolUnchanged++;
            } else if (rigid != null) {
                s.intra[m] = rigid.intra[m];
                numMolRigid++;
            } else if (partial != null) {
                double[] rOld = partial.r;
                int nMoved = 0;
                for (int i : cores) nMoved += markMoved(i, rOld, rNew, nMoved);
                for (int i : others) nMoved += markMoved(i, rOld, rNew, nMoved);
                // grid of the unmoved non-core atoms, which are at the same positions in rOld and rNew
                int nUnmoved = 0;
                for (int j : others) if (!moved[j]) unmovedList[nUnmoved++] = j;
                buildGrid(rNew, unmovedList, nUnmoved);
                double du = 0, scale = Math.abs(partial.intra[m]);
                for (int k = 0; k < nMoved; k++) {
                    int i = movedList[k];
                    du += uMovedLocal(rNew, rNew, i, cores, others, nMoved);
                    scale += absSum;
                    du -= uMovedLocal(rOld, rNew, i, cores, others, nMoved);
                    scale += absSum;
                }
                for (int k = 0; k < nMoved; k++) moved[movedList[k]] = false;
                double u = partial.intra[m] + du;
                // written so that NaN or infinity (e.g. infinite energies from coincident atoms,
                // giving inf - inf) also fails the test
                if (!(scale <= 1e4 * (1 + Math.abs(u)))) {
                    // large terms cancel (e.g. the stored or new configuration has overlapping
                    // atoms), so the update would lose precision: recompute instead
                    s.intra[m] = intraFull(rNew, m);
                    numMolFull++;
                } else {
                    s.intra[m] = u;
                    numMolPartial++;
                    sumMovedIntra += nMoved;
                }
            } else {
                s.intra[m] = intraFull(rNew, m);
                numMolFull++;
            }
        }
        computeInterIncremental(s);
    }

    private static int countMoved(double[] rOld, double[] rNew, int[] cores, int[] others) {
        int n = 0;
        for (int i : cores) if (atomMoved(rOld, rNew, i)) n++;
        for (int i : others) if (atomMoved(rOld, rNew, i)) n++;
        return n;
    }

    /**
     * True if atom i moved by more than round-off: a coordinate changed by more than
     * MOVE_TOLERANCE * (1 + |x|).  Undoing a rejected rotation (by rotating back) restores
     * positions only to round-off, and those atoms should count as unmoved.
     */
    private static boolean atomMoved(double[] rOld, double[] rNew, int i) {
        int i3 = 3 * i;
        for (int d = 0; d < 3; d++) {
            double a = rOld[i3 + d], b = rNew[i3 + d];
            if (a != b && Math.abs(a - b) > MOVE_TOLERANCE * (1 + Math.abs(a))) return true;
        }
        return false;
    }

    /** relative change in a coordinate below which an atom counts as unmoved */
    public static final double MOVE_TOLERANCE = 1e-12;

    /**
     * Returns true if molecule m moved rigidly from rOld to rNew: every atom's distances to four
     * non-coplanar reference atoms (chosen in rOld) are unchanged to within round-off.  Returns
     * false if the molecule is too close to planar for the test.
     */
    protected boolean isRigid(double[] rOld, double[] rNew, int m) {
        int[] all = new int[molCores[m].length + molOthers[m].length];
        System.arraycopy(molCores[m], 0, all, 0, molCores[m].length);
        System.arraycopy(molOthers[m], 0, all, molCores[m].length, molOthers[m].length);
        if (all.length < 4) return false;
        int r0 = all[0];
        // r1: farthest from r0
        int r1 = -1;
        double best = -1;
        for (int i : all) {
            double d = dist2(rOld, i, r0);
            if (d > best) {
                best = d;
                r1 = i;
            }
        }
        double size2 = best;
        if (size2 == 0) return false;
        // r2: farthest from the line r0-r1
        double[] e1 = unit(rOld, r0, r1);
        int r2 = -1;
        best = -1;
        for (int i : all) {
            double d = perp2(rOld, i, r0, e1);
            if (d > best) {
                best = d;
                r2 = i;
            }
        }
        if (best < 1e-6 * size2) return false;
        // r3: farthest from the plane r0-r1-r2
        double[] v2 = new double[3];
        for (int d = 0; d < 3; d++) v2[d] = rOld[3 * r2 + d] - rOld[3 * r0 + d];
        double[] nrm = {e1[1] * v2[2] - e1[2] * v2[1], e1[2] * v2[0] - e1[0] * v2[2], e1[0] * v2[1] - e1[1] * v2[0]};
        double nl = Math.sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
        int r3 = -1;
        best = -1;
        for (int i : all) {
            double h = 0;
            for (int d = 0; d < 3; d++) h += (rOld[3 * i + d] - rOld[3 * r0 + d]) * nrm[d];
            h = Math.abs(h) / nl;
            if (h > best) {
                best = h;
                r3 = i;
            }
        }
        if (best * best < 1e-6 * size2) return false;
        int[] refs = {r0, r1, r2, r3};
        double tol = 1e-10 * size2;
        for (int i : all) {
            for (int ref : refs) {
                if (Math.abs(dist2(rNew, i, ref) - dist2(rOld, i, ref)) > tol) return false;
            }
        }
        return true;
    }

    private static double dist2(double[] r, int i, int j) {
        double dx = r[3 * i] - r[3 * j], dy = r[3 * i + 1] - r[3 * j + 1], dz = r[3 * i + 2] - r[3 * j + 2];
        return dx * dx + dy * dy + dz * dz;
    }

    private static double[] unit(double[] r, int i, int j) {
        double[] e = {r[3 * j] - r[3 * i], r[3 * j + 1] - r[3 * i + 1], r[3 * j + 2] - r[3 * i + 2]};
        double l = Math.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        for (int d = 0; d < 3; d++) e[d] /= l;
        return e;
    }

    /** squared distance of atom i from the line through atom i0 along unit vector e */
    private static double perp2(double[] r, int i, int i0, double[] e) {
        double[] v = {r[3 * i] - r[3 * i0], r[3 * i + 1] - r[3 * i0 + 1], r[3 * i + 2] - r[3 * i0 + 2]};
        double p = v[0] * e[0] + v[1] * e[1] + v[2] * e[2];
        return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - p * p;
    }

    private int markMoved(int i, double[] rOld, double[] rNew, int nMoved) {
        if (atomMoved(rOld, rNew, i)) {
            moved[i] = true;
            movedList[nMoved] = i;
            return 1;
        }
        return 0;
    }

    /**
     * Energy of moved atom i with the atoms of its molecule (positions r), counting each
     * moved-moved pair once (from its lower-index atom).  Pairs involving a core and pairs of two
     * moved atoms are summed directly; pairs of a moved non-core atom with unmoved non-core atoms
     * come from the grid of unmoved atoms (built by the caller, positions rGrid; the unmoved atoms
     * are at the same positions in the old and new configurations, to round-off).  Sets absSum to the sum of the
     * absolute values of the terms.
     */
    private double uMovedLocal(double[] r, double[] rGrid, int i, int[] cores, int[] others, int nMoved) {
        double u = 0, a = 0;
        for (int j : cores) {
            if (j == i || (moved[j] && j < i)) continue;
            double v = uCore(r, i, j);
            u += v;
            a += Math.abs(v);
        }
        if (isCore[i]) {
            // a moved core: every non-core atom of the molecule, directly
            for (int j : others) {
                if (moved[j] && j < i) continue;
                double v = uCore(r, i, j);
                u += v;
                a += Math.abs(v);
            }
        } else {
            for (int k = 0; k < nMoved; k++) {
                int j = movedList[k];
                if (j <= i || isCore[j]) continue;
                double v = uOther(r, i, j);
                u += v;
                a += Math.abs(v);
            }
            u += gridEnergy(r[3 * i], r[3 * i + 1], r[3 * i + 2], atomType[i], rGrid);
            a += gridAbs;
        }
        absSum = a;
        return u;
    }

    /** sum of |u| over the terms of the last uMovedLocal call */
    private double absSum;

    /**
     * Interaction energies between all pairs of molecules, in full.
     */
    protected void computeInter(State s) {
        for (int m1 = 0; m1 < nMol; m1++) {
            for (int m2 = m1 + 1; m2 < nMol; m2++) {
                double u = interPair(s.r, m1, m2);
                s.inter[m1][m2] = s.inter[m2][m1] = u;
            }
        }
    }

    /**
     * Interaction energies for state s from the stored calculations: reused when neither molecule
     * moved; updated from the moved atoms when one molecule is unchanged and the other partly
     * moved; otherwise (or when the update would lose precision) recomputed in full.
     */
    protected void computeInterIncremental(State s) {
        double[] rNew = s.r;
        for (int m1 = 0; m1 < nMol; m1++) {
            for (int m2 = m1 + 1; m2 < nMol; m2++) {
                State unchanged = null, partial = null;
                int partialMol = -1, partialMoved = Integer.MAX_VALUE;
                for (State c : states) {
                    if (c == s || c.id < 0) continue;
                    int n1 = countMoved(c.r, rNew, molCores[m1], molOthers[m1]);
                    int n2 = countMoved(c.r, rNew, molCores[m2], molOthers[m2]);
                    if (n1 == 0 && n2 == 0) {
                        unchanged = c;
                        break;
                    }
                    int size1 = molCores[m1].length + molOthers[m1].length;
                    int size2 = molCores[m2].length + molOthers[m2].length;
                    int mMoved = -1, n = 0;
                    if (n2 == 0 && n1 < size1) {
                        mMoved = m1;
                        n = n1;
                    } else if (n1 == 0 && n2 < size2) {
                        mMoved = m2;
                        n = n2;
                    }
                    if (mMoved < 0) continue;
                    if (n < partialMoved || (n == partialMoved && Math.abs(c.inter[m1][m2]) < Math.abs(partial.inter[m1][m2]))) {
                        partial = c;
                        partialMol = mMoved;
                        partialMoved = n;
                    }
                }
                double u;
                if (unchanged != null) {
                    u = unchanged.inter[m1][m2];
                    numInterUnchanged++;
                } else if (partial != null) {
                    int ma = partialMol, mb = partialMol == m1 ? m2 : m1;
                    double[] rOld = partial.r;
                    int nMoved = 0;
                    for (int i : molCores[ma]) nMoved += markMoved(i, rOld, rNew, nMoved);
                    for (int i : molOthers[ma]) nMoved += markMoved(i, rOld, rNew, nMoved);
                    for (int k = 0; k < nMoved; k++) moved[movedList[k]] = false;
                    // molecule mb is unchanged: its non-core atoms in a grid, queried at the old
                    // and new positions of each moved atom of molecule ma
                    buildGrid(rNew, molOthers[mb], molOthers[mb].length);
                    double du = 0, scale = Math.abs(partial.inter[m1][m2]);
                    for (int k = 0; k < nMoved; k++) {
                        int i = movedList[k];
                        du += uWithMolecule(rNew, rNew, i, mb);
                        scale += absSum;
                        du -= uWithMolecule(rOld, rNew, i, mb);
                        scale += absSum;
                    }
                    u = partial.inter[m1][m2] + du;
                    // NaN or infinity (e.g. coincident atoms) also fails the test
                    if (!(scale <= 1e4 * (1 + Math.abs(u)))) {
                        // large terms cancel (overlapping molecules before or after): recompute
                        u = interPair(rNew, m1, m2);
                        numInterFull++;
                    } else {
                        numInterPartial++;
                        sumMovedInter += nMoved;
                    }
                } else {
                    u = interPair(rNew, m1, m2);
                    numInterFull++;
                }
                s.inter[m1][m2] = s.inter[m2][m1] = u;
            }
        }
    }

    /**
     * Energy of atom i (positions r) with every atom of molecule mb, whose non-core atoms are in
     * the grid (positions rGrid).  Sets absSum to the sum of the absolute values of the terms.
     */
    private double uWithMolecule(double[] r, double[] rGrid, int i, int mb) {
        double u = 0, a = 0;
        for (int j : molCores[mb]) {
            double v = uCore(r, i, j);
            u += v;
            a += Math.abs(v);
        }
        if (isCore[i]) {
            for (int j : molOthers[mb]) {
                double v = uCore(r, i, j);
                u += v;
                a += Math.abs(v);
            }
        } else {
            u += gridEnergy(r[3 * i], r[3 * i + 1], r[3 * i + 2], atomType[i], rGrid);
            a += gridAbs;
        }
        absSum = a;
        return u;
    }

    protected double interPair(double[] r, int m1, int m2) {
        double u = 0;
        // core atoms of each molecule with every atom of the other
        for (int i : molCores[m1]) {
            for (int j : molCores[m2]) u += uCore(r, i, j);
            for (int j : molOthers[m2]) u += uCore(r, i, j);
        }
        for (int i : molCores[m2]) {
            for (int j : molOthers[m1]) u += uCore(r, i, j);
        }
        // non-core pairs: grid of molecule m1's non-core atoms, queried by molecule m2's
        int[] a1 = molOthers[m1];
        buildGrid(r, a1, a1.length);
        for (int j : molOthers[m2]) {
            u += gridEnergy(r[3 * j], r[3 * j + 1], r[3 * j + 2], atomType[j], r);
        }
        return u;
    }

    /**
     * Puts atoms list[0..n-1] (non-core), at positions r, into a grid covering only their bounding
     * box, with cells at least rc wide (or a single cell in a direction narrower than rc).
     */
    protected void buildGrid(double[] r, int[] list, int n) {
        gridEmpty = n == 0;
        if (gridEmpty) return;
        for (int d = 0; d < 3; d++) {
            gMin[d] = Double.POSITIVE_INFINITY;
            gMax[d] = Double.NEGATIVE_INFINITY;
        }
        for (int k = 0; k < n; k++) {
            int i3 = 3 * list[k];
            for (int d = 0; d < 3; d++) {
                double x = r[i3 + d];
                if (x < gMin[d]) gMin[d] = x;
                if (x > gMax[d]) gMax[d] = x;
            }
        }
        for (int d = 0; d < 3; d++) {
            double extent = gMax[d] - gMin[d];
            gNc[d] = Math.max(1, Math.min(maxCellsPerDim, (int) (extent / rc)));
            gCellInv[d] = extent > 0 ? gNc[d] / extent : 0;
        }
        int totalCells = gNc[0] * gNc[1] * gNc[2];
        if (cellHead.length < totalCells) cellHead = new int[totalCells];
        java.util.Arrays.fill(cellHead, 0, totalCells, -1);
        for (int k = 0; k < n; k++) {
            int i = list[k];
            int c = cellIndex(r, i, gMin, gCellInv, gNc);
            cellNext[i] = cellHead[c];
            cellHead[c] = i;
        }
    }

    /**
     * Sum of the pair energies of a non-core atom of type ti at (x, y, z) with the atoms in the
     * grid (whose positions are in rGrid).  Also sets gridAbs to the sum of the absolute values.
     */
    protected double gridEnergy(double x, double y, double z, int ti, double[] rGrid) {
        gridAbs = 0;
        if (gridEmpty || x < gMin[0] - rc || x > gMax[0] + rc || y < gMin[1] - rc || y > gMax[1] + rc
                || z < gMin[2] - rc || z > gMax[2] + rc) return 0;
        int cx = clampCell((x - gMin[0]) * gCellInv[0], gNc[0]);
        int cy = clampCell((y - gMin[1]) * gCellInv[1], gNc[1]);
        int cz = clampCell((z - gMin[2]) * gCellInv[2], gNc[2]);
        IPotential2[] pi = atomPotentials[ti];
        double u = 0, a = 0;
        for (int ix = Math.max(0, cx - 1); ix <= Math.min(gNc[0] - 1, cx + 1); ix++) {
            for (int iy = Math.max(0, cy - 1); iy <= Math.min(gNc[1] - 1, cy + 1); iy++) {
                for (int iz = Math.max(0, cz - 1); iz <= Math.min(gNc[2] - 1, cz + 1); iz++) {
                    for (int j = cellHead[(iz * gNc[1] + iy) * gNc[0] + ix]; j > -1; j = cellNext[j]) {
                        int j3 = 3 * j;
                        double dx = rGrid[j3] - x, dy = rGrid[j3 + 1] - y, dz = rGrid[j3 + 2] - z;
                        double r2 = dx * dx + dy * dy + dz * dz;
                        if (r2 >= rc2) continue;
                        IPotential2 p = pi[atomType[j]];
                        if (p == null) continue;
                        double v = p.u(r2);
                        u += v;
                        a += Math.abs(v);
                    }
                }
            }
        }
        gridAbs = a;
        return u;
    }

    private static int clampCell(double x, int n) {
        int c = (int) x;
        return c < 0 ? 0 : (c >= n ? n - 1 : c);
    }

    private static int cellIndex(double[] r, int i, double[] min, double[] cellInv, int[] nc) {
        int cx = clampCell((r[3 * i] - min[0]) * cellInv[0], nc[0]);
        int cy = clampCell((r[3 * i + 1] - min[1]) * cellInv[1], nc[1]);
        int cz = clampCell((r[3 * i + 2] - min[2]) * cellInv[2], nc[2]);
        return (cz * nc[1] + cy) * nc[0] + cx;
    }

    /** pair energy when at least one atom is a core (no cutoff) */
    private double uCore(double[] r, int i, int j) {
        IPotential2 p = atomPotentials[atomType[i]][atomType[j]];
        if (p == null) return 0;
        int i3 = 3 * i, j3 = 3 * j;
        double dx = r[j3] - r[i3], dy = r[j3 + 1] - r[i3 + 1], dz = r[j3 + 2] - r[i3 + 2];
        return p.u(dx * dx + dy * dy + dz * dz);
    }

    /** pair energy of two non-core atoms (zero beyond the cutoff) */
    private double uOther(double[] r, int i, int j) {
        int i3 = 3 * i, j3 = 3 * j;
        double dx = r[j3] - r[i3], dy = r[j3 + 1] - r[i3 + 1], dz = r[j3 + 2] - r[i3 + 2];
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= rc2) return 0;
        IPotential2 p = atomPotentials[atomType[i]][atomType[j]];
        return p == null ? 0 : p.u(r2);
    }

    private void copyPositions(double[] r) {
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < nAtoms; i++) {
            Vector p = atoms.get(i).getPosition();
            r[3 * i] = p.getX(0);
            r[3 * i + 1] = p.getX(1);
            r[3 * i + 2] = p.getX(2);
        }
    }

    private void compare(State s, State check) {
        for (int m1 = 0; m1 < nMol; m1++) {
            if (differ(s.intra[m1], check.intra[m1])) {
                throw new RuntimeException("incremental intramolecular energy of molecule " + m1 + " is " + s.intra[m1]
                        + ", full recalculation gives " + check.intra[m1]);
            }
            for (int m2 = m1 + 1; m2 < nMol; m2++) {
                if (differ(s.inter[m1][m2], check.inter[m1][m2])) {
                    throw new RuntimeException("interaction energy of molecules " + m1 + " and " + m2 + " is "
                            + s.inter[m1][m2] + ", full recalculation gives " + check.inter[m1][m2]);
                }
            }
        }
    }

    private boolean differ(double a, double b) {
        return Math.abs(a - b) > debugTolerance * (1 + Math.abs(a) + Math.abs(b));
    }

    public PotentialCompute makePotentialComputeIntra() {
        return new PotentialComputeIntra();
    }

    /**
     * PotentialCompute giving the total intramolecular energy, for the moves that need it.
     */
    public class PotentialComputeIntra implements PotentialCompute {
        protected double lastEnergy;

        public void init() {
        }

        public Vector[] getForces() {
            return new Vector[0];
        }

        public double getLastVirial() {
            return 0;
        }

        public double getLastEnergy() {
            return lastEnergy;
        }

        public void updateAtom(IAtom atom) {
        }

        public double computeAll(boolean doForces, PotentialCallback pc) {
            double uintra = 0;
            for (IMolecule m : box.getMoleculeList()) {
                uintra += energy(m, m);
            }
            lastEnergy = uintra;
            return uintra;
        }

        public double computeOneOld(IAtom iAtom) {
            return 0;
        }

        public double computeOne(IAtom iAtom) {
            return 0;
        }

        public double computeManyAtomsOld(IAtom... atoms) {
            return 0;
        }

        public double computeManyAtoms(IAtom... atoms) {
            return 0;
        }

        public void processAtomU(double fac) {
        }

        public IntegratorListener makeIntegratorListener() {
            return new IntegratorListener() {
            };
        }
    }
}

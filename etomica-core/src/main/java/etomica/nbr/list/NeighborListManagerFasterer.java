package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.cell.NeighborCellManagerFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.Potential2Soft;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;

public class NeighborListManagerFasterer {
    private final NeighborCellManagerFasterer cellManager;
    private final Potential2Soft[][] pairPotentials;
    private final Box box;
    private final BondingInfo bondingInfo;
    private final boolean isPureAtoms;
    private final double[] maxR2, maxR2Unsafe;
    private final int numAtomTypes;
    private final Space space;
    public int[] numAtomNbrsUp, numAtomNbrsDn;
    // consider 1D array since Java sucks
    public int[][] nbrs;
    public Vector[][] nbrBoxOffsets;
    private double nbrRange;
    private double safetyFac = 0.45;
    private boolean onlyUpNbrs = true;
    private int maxNab;
    private Vector[] oldAtomPositions;

    public NeighborListManagerFasterer(Simulation sim, Box box, int cellRange, double nbrRange, Potential2Soft[][] pairPotentials, BondingInfo bondingInfo) {
        this.box = box;
        this.space = box.getSpace();
        this.bondingInfo = bondingInfo;
        this.cellManager = new NeighborCellManagerFasterer(box, cellRange);
        this.setNeighborRange(nbrRange);
        this.pairPotentials = pairPotentials;
        numAtomTypes = sim.getAtomTypeCount();
        maxR2 = new double[numAtomTypes];
        maxR2Unsafe = new double[numAtomTypes];
        oldAtomPositions = new Vector[0];
        nbrBoxOffsets = new Vector[0][0];
        nbrs = new int[0][0];
        numAtomNbrsUp = new int[0];
        isPureAtoms = sim.getSpeciesList().stream().allMatch(s -> s.getLeafAtomCount() == 1);
    }

    public double getNeighborRange() {
        return nbrRange;
    }

    public void setNeighborRange(double nbrRange) {
        this.nbrRange = nbrRange;
        cellManager.setPotentialRange(nbrRange);
    }

    public void setDoDownNeighbors(boolean doDown) {
        this.onlyUpNbrs = !doDown;
    }

    public void init() {
        cellManager.init();
        for (int i = 0; i < numAtomTypes; i++) {
            maxR2Unsafe[i] = maxR2[i] = 1e100;
            for (int j = 0; j < numAtomTypes; j++) {
                if (pairPotentials[i][j] == null) continue;
                double rc = pairPotentials[i][j].getRange();
                double maxDrUnsafe = (nbrRange - rc) * 0.5;
                double x = maxDrUnsafe * maxDrUnsafe;
                if (maxR2Unsafe[i] < x) continue;
                maxR2Unsafe[i] = x;
                double maxDr = maxDrUnsafe / 0.5 * safetyFac;
                maxR2[i] = maxDr * maxDr;
            }
        }
        reset();
    }

    public void reset() {
        IAtomList atoms = box.getLeafList();
        int boxNumAtoms = atoms.size();
        if (boxNumAtoms == 0) return;
        boolean moreAtoms = boxNumAtoms > oldAtomPositions.length;
        if (moreAtoms) {
            oldAtomPositions = new Vector[boxNumAtoms];
            for (int i = 0; i < boxNumAtoms; i++) oldAtomPositions[i] = box.getSpace().makeVector();
        }
        for (int i = 0; i < boxNumAtoms; i++) {
            Vector ri = atoms.get(i).getPosition();
            box.getBoundary().nearestImage(ri);
            oldAtomPositions[i].E(ri);
        }

        cellManager.assignCellAll();
        boolean forceReallocNbrs = false;
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();

        while (true) {
            if (moreAtoms) {
                numAtomNbrsUp = new int[boxNumAtoms];
                if (!onlyUpNbrs) numAtomNbrsDn = new int[boxNumAtoms];
            }
            // forceReallocNbrs can be used to force reallocation when max # of nbrs is too small
            if (moreAtoms || forceReallocNbrs) {
                maxNab *= 1.2;
                if (maxNab == 0) maxNab = 5;
                nbrs = new int[boxNumAtoms][maxNab];
                nbrBoxOffsets = new Vector[boxNumAtoms][maxNab];
                forceReallocNbrs = false;
            }

            for (int i = 0; i < boxNumAtoms; i++) numAtomNbrsUp[i] = 0;

            double rc2 = nbrRange * nbrRange;
            int tooMuch = 0;
            for (int i = 0; i < boxNumAtoms; i++) {
                IAtom iAtom = atoms.get(i);
                int j = i;
                int iCell = atomCell[i];
                Vector jbo = boxOffsets[iCell];
                Potential2Soft[] iPotentials = pairPotentials[iAtom.getType().getIndex()];
                while ((j = cellNextAtom[j]) > -1) {
                    IAtom jAtom = atoms.get(j);
                    tooMuch += checkNbrPair(i, j, iAtom, jAtom, rc2, jbo, iPotentials);
                }
                for (int cellOffset : cellOffsets) {
                    int jCell = iCell + cellOffset;
                    jbo = boxOffsets[jCell];
                    jCell = wrapMap[jCell];
                    for (j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                        IAtom jAtom = atoms.get(j);
                        tooMuch += checkNbrPair(i, j, iAtom, jAtom, rc2, jbo, iPotentials);
                    }
                }
                if (tooMuch > 0) {
                    if (Debug.ON) {
                        System.out.println("maxNab " + maxNab + " => " + (maxNab + tooMuch));
                    }
                    maxNab += tooMuch;
                    forceReallocNbrs = true;
                    break;
                }
            }
            if (tooMuch > 0) continue;
            if (!onlyUpNbrs) {
                for (int i = 0; i < boxNumAtoms; i++) numAtomNbrsDn[i] = 0;
                outerDn:
                for (int i = 0; i < boxNumAtoms; i++) {
                    int iNumNbrs = numAtomNbrsUp[i];
                    int[] iNbrs = nbrs[i];
                    for (int j = 0; j < iNumNbrs; j++) {
                        int jj = iNbrs[j];
                        if (numAtomNbrsDn[jj] + numAtomNbrsUp[jj] >= maxNab) {
                            if (Debug.ON) {
                                System.out.println("maxNab " + maxNab + " => " + (maxNab * 4 / 3));
                            }
                            maxNab = (maxNab * 4) / 3;
                            forceReallocNbrs = true;
                            break outerDn;
                        }
                        nbrs[jj][maxNab - 1 - numAtomNbrsDn[jj]] = i;
                        numAtomNbrsDn[jj]++;
                    }
                }
                if (forceReallocNbrs) continue;
            }
            break;
        }
    }

    private int checkNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, double rc2, Vector jbo, Potential2Soft[] iPotentials) {
        if (iPotentials[jAtom.getType().getIndex()] == null) return 0;

        if (bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom)) return 0;

        Vector dr = space.makeVector();
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double r2 = dr.squared();
        if (r2 > rc2) return 0;
        if (numAtomNbrsUp[i] >= maxNab) return 1;
        nbrs[i][numAtomNbrsUp[i]] = j;
        nbrBoxOffsets[i][numAtomNbrsUp[i]] = jbo;
        numAtomNbrsUp[i]++;
        return 0;
    }

    public void checkUpdateNbrs() {
        IAtomList atoms = box.getLeafList();
        int boxNumAtoms = atoms.size();
        boolean needsUpdate = false;
        for (int i = 0; i < boxNumAtoms; i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            double r2 = ri.Mv1Squared(oldAtomPositions[i]);
            if (r2 > maxR2[iAtom.getType().getIndex()]) {
                if (Debug.ON && safetyFac > 0 && r2 > maxR2Unsafe[iAtom.getType().getIndex()]) {
                    System.err.println(iAtom + " drifted into unsafe zone before nbr update");
                    needsUpdate = true;
                } else {
                    reset();
                    return;
                }
            }
        }
        if (Debug.ON && needsUpdate) reset();
    }
}

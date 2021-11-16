package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.BondingInfo;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2Soft;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.Debug;

import java.util.ArrayList;
import java.util.List;

public class NeighborListManager implements NeighborManager, NeighborManager.NeighborEventSource, IntegratorListener {
    private final NeighborCellManager cellManager;
    private Potential2Soft[][] pairPotentials;
    protected final Box box;
    protected final BondingInfo bondingInfo;
    protected final boolean isPureAtoms;
    private final double[] maxR2, maxR2Unsafe;
    private final int numAtomTypes;
    protected final Space space;
    private final NeighborIteratorList neighborIterator;
    public int[] numAtomNbrsUp, numAtomNbrsDn;
    // consider 1D array since Java sucks
    public int[][] nbrs;
    public Vector[][] nbrBoxOffsets;
    protected double nbrRange;
    private double safetyFac = 0.4;
    private boolean onlyUpNbrs = true;
    protected int maxNab;
    private Vector[] oldAtomPositions;
    private final List<INeighborListener> listeners;
    private int numUnsafe = -1;
    private double minR2;

    public NeighborListManager(SpeciesManager sm, Box box, int cellRange, double nbrRange, BondingInfo bondingInfo) {
        this.box = box;
        this.space = box.getSpace();
        this.bondingInfo = bondingInfo;
        this.cellManager = new NeighborCellManager(sm, box, cellRange, bondingInfo);
        this.setNeighborRange(nbrRange);
        numAtomTypes = sm.getAtomTypeCount();
        maxR2 = new double[numAtomTypes];
        maxR2Unsafe = new double[numAtomTypes];
        oldAtomPositions = new Vector[0];
        nbrBoxOffsets = new Vector[0][0];
        nbrs = new int[0][0];
        numAtomNbrsUp = new int[0];
        isPureAtoms = sm.isPureAtoms();
        this.neighborIterator = new NeighborIteratorList(this, box);
        listeners = new ArrayList<>();
    }

    public NeighborCellManager getCellManager() {
        return cellManager;
    }

    public Box getBox() {
        return box;
    }

    @Override
    public BondingInfo getBondingInfo() {
        return bondingInfo;
    }

    public void addListener(INeighborListener newListener) {
        listeners.add(newListener);
    }

    @Override
    public void setPairPotentials(Potential2Soft[][] potentials) {
        this.pairPotentials = potentials;
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

    @Override
    public NeighborIterator makeNeighborIterator() {
        return this.neighborIterator;
    }

    public void setSafetyFac(double newSafetyFac) {
        safetyFac = newSafetyFac;
        init();
    }

    public void init() {
        if (pairPotentials == null) return;
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

    @Override
    public IntegratorListener makeIntegratorListener() {
        return this;
    }

    @Override
    public void integratorInitialized(IntegratorEvent e) {
        init();
    }

    @Override
    public void integratorStepStarted(IntegratorEvent e) {
        checkUpdateNbrs();
    }


    @Override
    public void updateAtom(IAtom atom) {
        this.cellManager.updateAtom(atom);
    }

    protected void realloc() {
        int boxNumAtoms = box.getLeafList().size();
        nbrs = new int[boxNumAtoms][maxNab];
        nbrBoxOffsets = new Vector[boxNumAtoms][maxNab];
    }

    public void reset() {
        Vector boxSize = box.getBoundary().getBoxSize();
        double minL = Double.POSITIVE_INFINITY;
        for (int i=0; i<boxSize.getD(); i++) {
            minL = Math.min(minL, boxSize.getX(i));
        }
        minR2 = (0.5*minL)*(0.5*minL);

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
            ri.PE(box.getBoundary().centralImage(ri));
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

        if (!moreAtoms && boxNumAtoms > nbrs.length) realloc();

        while (true) {
            if (moreAtoms) {
                numAtomNbrsUp = new int[boxNumAtoms];
                if (!onlyUpNbrs) numAtomNbrsDn = new int[boxNumAtoms];
            }
            // forceReallocNbrs can be used to force reallocation when max # of nbrs is too small
            if (moreAtoms || forceReallocNbrs) {
                if (forceReallocNbrs) {
                    maxNab *= 1.2;
                }
                if (maxNab == 0) maxNab = 5;
                realloc();
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
                IPotentialAtomic[] iPotentials = pairPotentials[iAtom.getType().getIndex()];
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
                        newDownNeighbor(jj, i, j, maxNab - 1 - numAtomNbrsDn[jj]);
                        // these offsets are actually opposite... we don't have a clean
                        // way to find the index of the opposite offset.
                        nbrBoxOffsets[jj][maxNab - 1 - numAtomNbrsDn[jj]] = nbrBoxOffsets[i][j];

                        numAtomNbrsDn[jj]++;
                    }
                }
                if (forceReallocNbrs) continue;
            }
            break;
        }
    }

    protected void newDownNeighbor(int j, int i, int upSlot, int downSlot) {
        nbrs[j][downSlot] = i;
    }

    protected int checkNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, double rc2, Vector jbo, IPotentialAtomic[] iPotentials) {
        if (iPotentials[jAtom.getType().getIndex()] == null) return 0;

        boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom);

        Vector dr = space.makeVector();
        Vector ri = iAtom.getPosition();
        Vector rj = jAtom.getPosition();
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double r2 = dr.squared();
        if (r2 > rc2 || (skipIntra && r2 < minR2)) return 0;
        return addAsNbrPair(i, j, iAtom, jAtom, jbo, iPotentials, dr);
    }

    protected int addAsNbrPair(int i, int j, IAtom iAtom, IAtom jAtom, Vector jbo, IPotentialAtomic[] iPotentials, Vector dr) {
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
        double thisMaxR2 = 0;
        boolean unsafe = false;
        for (int i = 0; i < boxNumAtoms; i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            double r2 = ri.Mv1Squared(oldAtomPositions[i]);
            thisMaxR2 = Math.max(r2, thisMaxR2);
            if (r2 > maxR2[iAtom.getType().getIndex()]) {
                if (safetyFac > 0 && r2 > maxR2Unsafe[iAtom.getType().getIndex()]) {
                    unsafe = true;
                    if (numUnsafe == -1) {
                        numUnsafe++;
                        System.out.println(iAtom + " drifted into unsafe zone before nbr update " + Math.sqrt(r2) + " > " + Math.sqrt(maxR2Unsafe[iAtom.getType().getIndex()]));
                    }
                }
                needsUpdate = true;
            }
        }

        if (needsUpdate) {
            if (unsafe) {
                numUnsafe++;
                if (numUnsafe == 1 || (Long.toString(numUnsafe).matches("10*"))) {
                    System.err.print("Atoms exceeded the safe neighbor limit");
                    if (numUnsafe > 1) {
                        System.err.print(" (" + numUnsafe + " times)");
                    }
                    System.err.println();
                }
            }
            reset();
            fireNeighborUpdateEvent();
        }
    }

    private void fireNeighborUpdateEvent() {
        for (INeighborListener l : listeners) {
            l.neighborListNeighborsUpdated();
        }
    }
}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.integrator;

import etomica.atom.*;
import etomica.box.Box;
import etomica.box.BoxMoleculeEvent;
import etomica.molecule.IMolecule;
import etomica.nbr.list.INeighborListener;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.IPotential1;
import etomica.potential.IPotential2;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.*;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.collections.Int2IntHash;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IntegratorHard extends IntegratorMD implements INeighborListener {

    protected final NeighborManagerHard neighborManager;
    protected final PotentialMasterBonding.FullBondingInfo bondingInfo;
    protected double[] collisionTimes, nullCollisionTimes;
    protected int[] collisionPartners;
    protected IPotential2[] collisionPotentials;
    protected Vector[] collisionVector;
    protected int[] eventBinsFirstAtom, eventBinsNextAtom, eventBinsPrevAtom;
    protected int[] collisionOldState;
    protected final IPotential2[][] pairPotentials;
    protected final IPotential1[] fieldPotentials;
    protected final double[] maxSigma;
    protected double tBase, tMax;
    protected int collisionUpdateInterval, collisionUpdateCountdown;
    protected long collisionCount, intervalCollisionCount;
    protected final NeighborIteratorHard neighborIterator;
    protected final List<CollisionListener> collisionListeners;
    private int lastBin = 0;
    private Int2IntHash[] bondState;
    private int[] fieldState;
    private final IAtom[] nextCollidingPair = new IAtom[2];
    private int nextColliderIndex;

    private long tDelete, tAdd, tNext, tSteps, tUpdate, tUp, tDown, tBump, tData, tCollect, tCollision, tAdvance, tNotStep;
    private final boolean writeTiming = false, verbose = false;
    private long lastTimingPrint = 0;
    private int handlingEvent = 0;

    private long nanoTime() {
        return writeTiming ? System.nanoTime() : 0;
    }

    private static PotentialCompute makeTotalCompute(Box box, IPotential2[][] pairPotentials, NeighborManagerHard neighborManager, IPotential1[] fieldPotentials, SpeciesManager sm, PotentialMasterBonding.FullBondingInfo bondingInfo) {
        int n = pairPotentials.length;
        IPotential2[][] softPotentials = new IPotential2[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                softPotentials[i][j] = (IPotential2) pairPotentials[i][j];
            }
        }
        List<PotentialCompute> computes = new ArrayList<>();
        PotentialComputePairGeneral pcPair = new PotentialComputePairGeneral(sm, softPotentials, box, neighborManager);
        computes.add(pcPair);
        if (fieldPotentials != null) {
            PotentialComputeField pcField = new PotentialComputeField(fieldPotentials, box);
            computes.add(pcField);
        }
        if (bondingInfo != null) {
            PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box, bondingInfo);
            computes.add(pmBonding);
        }
        if (computes.size() == 1) return pcPair;
        PotentialCompute[] computesArray = computes.toArray(new PotentialCompute[0]);
        return new PotentialComputeAggregate(computesArray);
    }

    /**
     * Constructs integrator with a default for non-isothermal sampling.
     *
     * @param random           random number generator used for initial velocities and some thermostats
     * @param timeStep         time step for integration
     * @param temperature      used by thermostat and/or to initialize velocities
     * @param box
     */
    public IntegratorHard(IPotential2[][] pairPotentials, NeighborManagerHard neighborManager, IRandom random, double timeStep, double temperature, Box box, SpeciesManager sm) {
        this(pairPotentials, neighborManager, random, timeStep, temperature, box, sm, null);
    }

    public IntegratorHard(IPotential2[][] pairPotentials, NeighborManagerHard neighborManager, IRandom random, double timeStep, double temperature, Box box, SpeciesManager sm, PotentialMasterBonding.FullBondingInfo bondingInfo) {
        this(pairPotentials, null, neighborManager, random, timeStep, temperature, box, sm, bondingInfo);
    }

    public IntegratorHard(IPotential2[][] pairPotentials, IPotential1[] fieldPotentials, NeighborManagerHard neighborManager, IRandom random, double timeStep, double temperature, Box box, SpeciesManager sm, PotentialMasterBonding.FullBondingInfo bondingInfo) {
        this(makeTotalCompute(box, pairPotentials, neighborManager, fieldPotentials, sm, bondingInfo), pairPotentials, fieldPotentials, neighborManager, random, timeStep, temperature, box, bondingInfo);
    }

    public IntegratorHard(PotentialCompute pc, IPotential2[][] pairPotentials, IPotential1[] fieldPotentials, NeighborManagerHard neighborManager, IRandom random, double timeStep, double temperature, Box box, PotentialMasterBonding.FullBondingInfo bondingInfo) {
        super(pc, random, timeStep, temperature, box);
        this.pairPotentials = pairPotentials;
        this.fieldPotentials = fieldPotentials;
        this.neighborManager = neighborManager;
        this.bondingInfo = bondingInfo;
        eventBinsFirstAtom = new int[10];
        collisionVector = new Vector[0];
        maxSigma = new double[pairPotentials.length];
        resizeArrays();
        // we want to store enough collisions that we don't have to recompute them
        // before the neighbor update
        collisionUpdateInterval = 5;
        tMax = collisionUpdateInterval * timeStep;
        neighborIterator = neighborManager.makeNeighborIterator();
        if (neighborManager instanceof NeighborManager.NeighborEventSource) {
            ((NeighborManager.NeighborEventSource) neighborManager).addListener(this);
        } else {
            for (int i = 0; i < box.getSpace().D(); i++) {
                if (box.getBoundary().getPeriodicity(i)) {
                    nullCollisionTimes = new double[box.getLeafList().size()];
                    break;
                }
            }
        }
        collisionListeners = new ArrayList<>(0);
    }

    private static PotentialCompute makeTotalCompute(PotentialComputePair pcPair, PotentialComputeField pcField) {
        if (pcPair == null) return pcField;
        if (pcField == null) return pcPair;
        return new PotentialComputeAggregate(pcField, pcPair);
    }

    public void setPairPotential(AtomType type1, AtomType type2, IPotential2 p2) {
        pairPotentials[type1.getIndex()][type2.getIndex()] = p2;
        pairPotentials[type2.getIndex()][type1.getIndex()] = p2;
        if (potentialCompute instanceof PotentialComputeAggregate) {
            for (PotentialCompute pc : ((PotentialComputeAggregate) potentialCompute).getPotentialComputes()) {
                if (pc instanceof PotentialComputePair) {
                    ((PotentialComputePair) pc).setPairPotential(type1, type2, p2);
                } else if (pc instanceof PotentialComputePairGeneral) {
                    ((PotentialComputePairGeneral) pc).setPairPotential(type1, type2, p2);
                }
            }
        } else if (potentialCompute instanceof PotentialComputePair) {
            ((PotentialComputePair) potentialCompute).setPairPotential(type1, type2, p2);
        } else if (potentialCompute instanceof PotentialComputePairGeneral) {
            ((PotentialComputePairGeneral) potentialCompute).setPairPotential(type1, type2, p2);
        }
    }

    public void setFieldPotential(AtomType type, IPotential1 pField) {
        fieldPotentials[type.getIndex()] = pField;
        if (potentialCompute instanceof PotentialComputeAggregate) {
            for (PotentialCompute pc : ((PotentialComputeAggregate) potentialCompute).getPotentialComputes()) {
                if (pc instanceof PotentialComputeField) {
                    ((PotentialComputeField) pc).setFieldPotential(type, pField);
                }
            }
        } else if (potentialCompute instanceof PotentialComputeField) {
            ((PotentialComputeField) potentialCompute).setFieldPotential(type, pField);
        }

    }

    public void setCollisionUpdateInterval(int newInterval) {
        collisionUpdateInterval = newInterval;
        tMax = collisionUpdateInterval * timeStep;
        if (initialized) reset();
    }

    public void setMaxCollisionDiameter(AtomType atomType, double maxSigma) {
        this.maxSigma[atomType.getIndex()] = maxSigma;
        if (nullCollisionTimes == null) {
            nullCollisionTimes = new double[box.getLeafList().size()];
        }
    }

    public void addCollisionListener(CollisionListener newListener) {
        collisionListeners.add(newListener);
    }

    public void removeCollisionListener(CollisionListener oldListener) {
        collisionListeners.remove(oldListener);
    }

    public void setTimeStep(double dt) {
        super.setTimeStep(dt);
        tMax = collisionUpdateInterval * timeStep;
        if (initialized) reset();
    }

    public void resetStepCount() {
        super.resetStepCount();
        lastTimingPrint = 0;
    }

    protected void resizeArrays() {
        int nAtoms = box.getLeafList().size();
        if (collisionVector.length < nAtoms) {
            collisionTimes = new double[nAtoms];
            collisionPartners = new int[nAtoms];
            eventBinsNextAtom = new int[nAtoms];
            eventBinsPrevAtom = new int[nAtoms];
            collisionOldState = new int[nAtoms];
            if (bondingInfo != null) {
                bondState = new Int2IntHash[nAtoms];
                for (int i = 0; i < nAtoms; i++) {
                    bondState[i] = new Int2IntHash(2);
                }
            }
            collisionPotentials = new IPotential2[nAtoms];
            collisionVector = new Vector[nAtoms];
            for (int i = 0; i < nAtoms; i++) {
                collisionVector[i] = space.makeVector();
            }
            if (fieldPotentials != null) {
                fieldState = new int[nAtoms];
            }
            if (nullCollisionTimes != null) {
                nullCollisionTimes = new double[nAtoms];
            }
        }
    }

    public void reset() {
        if (verbose) System.out.println("resetting integrator");
        super.reset();
        resizeArrays();
        neighborManager.init();

        computeAllCollisions();
    }

    protected long nPerBinSum = 0, nComputeAll = 0, nPerBinMax = 0;

    protected void computeAllCollisions() {
        if (!initialized) {
            // we got triggered because something major happened, but we're still getting set up
            // we'll get called from reset later
            return;
        }

        tBase = 0;
        IAtomList atoms = box.getLeafList();
        for (IAtom atom : atoms) {
            IAtomKinetic iAtom = (IAtomKinetic) atom;
            computeNullCollisionTime(iAtom, 0);
            int i = iAtom.getLeafIndex();
            if (bondingInfo != null) {

                IMolecule parentMolecule = iAtom.getParentGroup();
                Vector ri = iAtom.getPosition();
                bondingInfo.bondedPairPartners[iAtom.getParentGroup().getType().getIndex()].forEach((potential, partners) -> {
                    for (int[] pairIdx : partners[iAtom.getIndex()]) {
                        IAtomKinetic atom1 = (IAtomKinetic) parentMolecule.getChildList().get(pairIdx[0]);
                        IAtomKinetic atom2 = (IAtomKinetic) parentMolecule.getChildList().get(pairIdx[1]);
                        IAtomKinetic jAtom = atom1 == iAtom ? atom2 : atom1;
                        int j = jAtom.getLeafIndex();
                        if (j < i) continue; // no need to handle each pair twice
                        Vector rj = jAtom.getPosition();
                        Vector dr = space.makeVector();
                        dr.Ev1Mv2(rj, ri);
                        box.getBoundary().nearestImage(dr);
                        int s = potential.getState(iAtom, jAtom, dr);
                        bondState[i].put(j, s);
                        bondState[j].put(i, s);
                    }
                });
            }
            if (fieldPotentials != null && fieldPotentials[iAtom.getType().getIndex()] != null) {
                fieldState[i] = fieldPotentials[iAtom.getType().getIndex()].getState(iAtom);
            }
        }

        lastBin = 0;
        Arrays.fill(eventBinsFirstAtom, -1);
        Arrays.fill(collisionTimes, Double.POSITIVE_INFINITY);
        for (IAtom atom : box.getLeafList()) {
            collisionTimeUp((IAtomKinetic) atom, 0);
        }
        collisionUpdateCountdown = collisionUpdateInterval;
        int n = 0;
        int i = eventBinsFirstAtom[0];
        if (i >= 0) {
            n++;
            while ((i = eventBinsNextAtom[i]) > -1) {
                n++;
            }
        }
        nPerBinSum += n;
        nComputeAll++;
        nPerBinMax = Math.max(nPerBinMax, n);

        if (n > 50) {
            eventBinsFirstAtom = new int[((eventBinsFirstAtom.length + 1) * 3) / 2];
            if (verbose)
                System.out.println("bins increased to " + eventBinsFirstAtom.length + " -- " + (tMax / eventBinsFirstAtom.length) + " time per bin");
            Arrays.fill(eventBinsFirstAtom, -1);
            for (int j = 0; j < atoms.size(); j++) {
                if (collisionTimes[j] >= tMax) continue;
                int b = (int) (collisionTimes[j] / tMax * eventBinsFirstAtom.length);

                int oldFirst = eventBinsFirstAtom[b];
                eventBinsNextAtom[j] = oldFirst;
                eventBinsPrevAtom[j] = -1 - b;
                if (oldFirst >= 0) eventBinsPrevAtom[oldFirst] = j;
                eventBinsFirstAtom[b] = j;
            }

        }
    }

    protected int nextCollider() {
        long t1 = nanoTime();
        for (int jj = lastBin; jj < eventBinsFirstAtom.length; jj++) {
            int first = eventBinsFirstAtom[jj];
            if (first >= 0) {
                int i = first;
                double tFirst = collisionTimes[i];
                while ((i = eventBinsNextAtom[i]) > -1) {
                    if (collisionTimes[i] < tFirst) {
                        first = i;
                        tFirst = collisionTimes[i];
                    }
                }
                lastBin = jj;
                tNext += nanoTime() - t1;
                return first;
            }
        }
        lastBin = eventBinsFirstAtom.length;
        tNext += nanoTime() - t1;
        return -1;
    }

    protected void collectDownColliders(IAtom atom, AtomArrayList downColliders) {
        int i = atom.getLeafIndex();
        neighborIterator.iterDownNeighbors(atom.getLeafIndex(), (jAtom, rij, n) -> {
            if (collisionPartners[jAtom.getLeafIndex()] == i) {
                downColliders.add(jAtom);
            }
        });
        if (bondingInfo != null) {
            IMolecule parentMolecule = atom.getParentGroup();
            bondingInfo.bondedPairPartners[parentMolecule.getType().getIndex()].forEach((potential, partners) -> {
                for (int[] pairIdx : partners[atom.getIndex()]) {
                    IAtom atom1 = parentMolecule.getChildList().get(pairIdx[0]);
                    IAtom atom2 = parentMolecule.getChildList().get(pairIdx[1]);
                    IAtom jAtom = atom1 == atom ? atom2 : atom1;
                    // we get all bonded pairs containing i.  skip the "up" neighbors
                    if (jAtom.getLeafIndex() < atom.getLeafIndex() && collisionPartners[jAtom.getLeafIndex()] == i) {
                        downColliders.add(jAtom);
                    }
                }
            });
        }

    }

    protected void removeCollision(int i) {
        long t1 = nanoTime();
        int prev = eventBinsPrevAtom[i];
        int next = eventBinsNextAtom[i];
        if (prev < 0) {
            eventBinsFirstAtom[-prev - 1] = next;
        } else {
            eventBinsNextAtom[prev] = next;
        }
        if (next >= 0) {
            eventBinsPrevAtom[next] = prev;
        }
        tDelete += nanoTime() - t1;
    }

    protected void computeNullCollisionTime(IAtomKinetic atom, double falseTime) {
        if (nullCollisionTimes == null) return;
        Vector v = atom.getVelocity();
        Vector dim = box.getBoundary().getBoxSize();
        double tmin = Double.POSITIVE_INFINITY;
        double d2 = 2.0 * maxSigma[atom.getType().getIndex()];
        int D = dim.getD();
        // 4*(L/4 - sigma/2) = L - 2sigma
        for (int i = 0; i < D; i++) {
            double t = (dim.getX(i) - d2) / v.getX(i);
            t = Math.abs(t);
            tmin = Math.min(tmin, t);
        }
        nullCollisionTimes[atom.getLeafIndex()] = tBase + falseTime + 0.25 * tmin;
    }

    protected void collisionTimeDown(IAtomKinetic iAtom, double falseTime) {
        long t1cd = nanoTime();
        Vector iVelocity = iAtom.getVelocity();
        int iType = iAtom.getType().getIndex();
        int i = iAtom.getLeafIndex();
        MyConsumer myDown = (jAtom, rij, state, potential) -> {
            int j = jAtom.getLeafIndex();
            IAtomKinetic jjAtom = (IAtomKinetic) jAtom;
            Vector dv = space.makeVector();
            dv.Ev1Mv2(jjAtom.getVelocity(), iVelocity);
            if (neighborManager instanceof NeighborListManager) {
                rij.PEa1Tv1(falseTime, dv);
            }
            long t1 = nanoTime();
            double time = potential.collisionTime(iAtom, jjAtom, rij, dv, state, falseTime);
            tCollision += nanoTime() - t1;
            if (time < 0) {
                System.out.println("computed col down " + i + " with " + j + " at " + time);
                Vector dr = space.makeVector();
                dr.Ev1Mv2(jAtom.getPosition(), iAtom.getPosition());
                System.out.println("dr " + rij + " " + Math.sqrt(rij.squared()) + " bij " + rij.dot(dv) + " state " + state);
                System.out.println(j + " thought it would collide at " + collisionTimes[j] + " with " + collisionPartners[j]);
                if (nullCollisionTimes != null) {
                    System.out.println(i + " null time: " + nullCollisionTimes[i]);
                    System.out.println(j + " null time: " + nullCollisionTimes[j]);
                }
                potential.collisionTime(iAtom, jjAtom, rij, dv, state, falseTime);
                throw new RuntimeException("negative!");
            }

            double tcol = time + falseTime + tBase;

            if (tcol < Math.min(collisionTimes[j], tMax)) {
                if (collisionTimes[j] < tMax) {
                    // need to remove j from event bins
                    removeCollision(j);
                }
                collisionTimes[j] = tcol;
                collisionPartners[j] = i;
                collisionPotentials[j] = potential;
                collisionOldState[j] = state;
                rij.PEa1Tv1(time, dv);
                collisionVector[j].Ea1Tv1(-1, rij);

                long t1a = nanoTime();
                int b = (int) (collisionTimes[j] / tMax * eventBinsFirstAtom.length);

                int oldFirst = eventBinsFirstAtom[b];
                eventBinsNextAtom[j] = oldFirst;
                eventBinsPrevAtom[j] = -1 - b;
                if (oldFirst >= 0) eventBinsPrevAtom[oldFirst] = j;
                eventBinsFirstAtom[b] = j;
                tAdd += nanoTime() - t1a;
            }
        };
        NeighborConsumerHard ncDown = new NeighborConsumerHard() {
            @Override
            public void acceptHard(IAtom jAtom, Vector rij, int state) {
                myDown.accept(jAtom, rij, state, pairPotentials[iType][jAtom.getType().getIndex()]);
            }

            @Override
            public void accept(IAtom jAtom, Vector rij, int n) {

            }
        };
        neighborIterator.iterDownNeighbors(iAtom.getLeafIndex(), ncDown, falseTime);

        if (bondingInfo != null) {
            IMolecule parentMolecule = iAtom.getParentGroup();
            Vector ri = iAtom.getPosition();
            bondingInfo.bondedPairPartners[iAtom.getParentGroup().getType().getIndex()].forEach((potential, partners) -> {
                for (int[] pairIdx : partners[iAtom.getIndex()]) {
                    IAtom atom1 = parentMolecule.getChildList().get(pairIdx[0]);
                    IAtom atom2 = parentMolecule.getChildList().get(pairIdx[1]);
                    IAtom jAtom = atom1 == iAtom ? atom2 : atom1;
                    // we get all bonded pairs containing i.  skip the "up" neighbors
                    int j = jAtom.getLeafIndex();
                    if (j > i) continue;
                    Vector rj = jAtom.getPosition();
                    Vector dr = space.makeVector();
                    dr.Ev1Mv2(rj, ri);
                    box.getBoundary().nearestImage(dr);
                    myDown.accept(jAtom, dr, bondState[j].get(i), potential);
                }
            });
        }

        tDown += nanoTime() - t1cd;
    }

    protected void collisionTimeUp(IAtomKinetic iAtom, double falseTime) {
        long t1cu = nanoTime();
        int i = iAtom.getLeafIndex();
        if (collisionTimes[i] < tMax) {
            // need to remove i from event bins
            removeCollision(i);
        }
        collisionTimes[i] = Double.POSITIVE_INFINITY;
        Vector iVelocity = iAtom.getVelocity();
        int iType = iAtom.getType().getIndex();
        double[] tCheck = new double[]{tMax - tBase - falseTime};
        MyConsumer myUp = (jAtom, rij, state, potential) -> {
            IAtomKinetic jjAtom = (IAtomKinetic) jAtom;
            Vector dv = space.makeVector();
            dv.Ev1Mv2(jjAtom.getVelocity(), iVelocity);
            if (neighborManager instanceof NeighborListManager) {
                rij.PEa1Tv1(falseTime, dv);
            }
            long t1ct = nanoTime();
            int j = jAtom.getLeafIndex();
            double time = potential.collisionTime(iAtom, jjAtom, rij, dv, state, falseTime);
            tCollision += nanoTime() - t1ct;
            if (time < 0) {
                System.out.println("computed col up " + i + " with " + j + " at " + time);
                Vector r1 = box.getSpace().makeVector();
                r1.E(iAtom.getPosition());
                r1.PEa1Tv1(falseTime, iAtom.getVelocity());
                System.out.println(iAtom + " " + iAtom.getPosition() + " " + r1);
                r1.E(jAtom.getPosition());
                r1.PEa1Tv1(falseTime, jjAtom.getVelocity());
                System.out.println(jAtom + " " + jAtom.getPosition() + " " + r1);
                System.out.println("dr " + rij + "  |r| " + Math.sqrt(rij.squared()));
                System.out.println("dv " + dv + "  bij " + rij.dot(dv));
                System.out.println("state " + state);
                potential.collisionTime(iAtom, jjAtom, rij, dv, state, falseTime);
                throw new RuntimeException("Negative up collision time between " + i + " and " + j + " at " + time);
            }
            if (time < tCheck[0]) {
//                long t1 = nanoTime();
                collisionTimes[i] = time;
                tCheck[0] = time;
                collisionPartners[i] = j;
                collisionPotentials[i] = potential;
                rij.PEa1Tv1(time, dv);
                collisionVector[i].E(rij);
                collisionOldState[i] = state;
//                tAdd += nanoTime() - t1;
            }
        };

        NeighborConsumerHard up = new NeighborConsumerHard() {
            @Override
            public void acceptHard(IAtom jAtom, Vector rij, int state) {
                myUp.accept(jAtom, rij, state, pairPotentials[iType][jAtom.getType().getIndex()]);

            }

            @Override
            public void accept(IAtom jAtom, Vector rij, int n) {

            }
        };
        neighborIterator.iterUpNeighbors(i, up, falseTime);

        if (bondingInfo != null) {
            IMolecule parentMolecule = iAtom.getParentGroup();
            Vector ri = iAtom.getPosition();
            bondingInfo.bondedPairPartners[iAtom.getParentGroup().getType().getIndex()].forEach((potential, partners) -> {
                for (int[] pairIdx : partners[iAtom.getIndex()]) {
                    IAtom atom1 = parentMolecule.getChildList().get(pairIdx[0]);
                    IAtom atom2 = parentMolecule.getChildList().get(pairIdx[1]);
                    IAtom jAtom = atom1 == iAtom ? atom2 : atom1;
                    // we get all bonded pairs containing i.  skip the "down" neighbors
                    int j = jAtom.getLeafIndex();
                    if (j < i) continue;
                    Vector rj = jAtom.getPosition();
                    Vector dr = space.makeVector();
                    dr.Ev1Mv2(rj, ri);
                    box.getBoundary().nearestImage(dr);
                    myUp.accept(jAtom, dr, bondState[i].get(j), potential);
                }
            });
        }

        IPotential1 pField = fieldPotentials != null ? fieldPotentials[iType] : null;
        if (pField != null) {
            Vector rfalse = iAtom.getPosition();
            Vector r = Vector.d(rfalse.getD());
            r.E(rfalse);
            Vector v = iAtom.getVelocity();
            r.PEa1Tv1(falseTime, v);
            long t1ct = nanoTime();
            double time = pField.collisionTime(iAtom, r, v, fieldState[i], falseTime);
            tCollision += nanoTime() - t1ct;
            if (time < 0) {
                throw new RuntimeException("Negative up collision time between " + i + " and wall at " + time);
            }
            if (time < tCheck[0]) {
//                long t1 = nanoTime();
                collisionTimes[i] = time;
                tCheck[0] = time;
                collisionPartners[i] = -1; // partner = -1 means field
                collisionPotentials[i] = null;
                r.PEa1Tv1(time, v);
                collisionVector[i].E(r);
                collisionOldState[i] = fieldState[i];
//                tAdd += nanoTime() - t1;
            }
        }

        collisionTimes[i] += tBase + falseTime;

        if (nullCollisionTimes != null && nullCollisionTimes[i] < collisionTimes[i]) {
            // revert to null collision
            collisionTimes[i] = nullCollisionTimes[i];
            collisionPartners[i] = -2;
            collisionPotentials[i] = null;
        }

        long t1 = nanoTime();
        int b = (int) (collisionTimes[i] / tMax * eventBinsFirstAtom.length);

        if (b < eventBinsFirstAtom.length) {
            // collision time is less than max.  a keeper!
//            eventBinsNextAtom[i] = eventBinsFirstAtom[b];
//            eventBinsFirstAtom[b] = i;

            int oldFirst = eventBinsFirstAtom[b];
            eventBinsNextAtom[i] = oldFirst;
            eventBinsPrevAtom[i] = -1 - b;
            if (oldFirst >= 0) eventBinsPrevAtom[oldFirst] = i;
            eventBinsFirstAtom[b] = i;

        }
        tAdd += nanoTime() - t1;
        tUp += nanoTime() - t1cu;
    }

    protected void updatePair(IAtomKinetic atom1, IAtomKinetic atom2, double falseTime) {
        long t1 = nanoTime();
        AtomArrayList downColliders = new AtomArrayList(0);
        computeNullCollisionTime(atom1, falseTime);
        computeNullCollisionTime(atom2, falseTime);
        collectDownColliders(atom1, downColliders);
        collectDownColliders(atom2, downColliders);
        tCollect += nanoTime() - t1;
        for (IAtom iAtom : downColliders) {
            collisionTimeUp((IAtomKinetic) iAtom, falseTime);
        }

        // don't need to do up for atom1 because it was included above!
        collisionTimeDown(atom1, falseTime);

        collisionTimeUp(atom2, falseTime);
        collisionTimeDown(atom2, falseTime);
        tUpdate += nanoTime() - t1;
    }

    protected void updateAtom(IAtomKinetic atom, double falseTime) {
        long t1 = nanoTime();
        computeNullCollisionTime(atom, falseTime);
        AtomArrayList downColliders = new AtomArrayList(0);
        collectDownColliders(atom, downColliders);
        tCollect += nanoTime() - t1;
        for (IAtom iAtom : downColliders) {
            collisionTimeUp((IAtomKinetic) iAtom, falseTime);
        }

        collisionTimeDown(atom, falseTime);
        collisionTimeUp(atom, falseTime);

        tUpdate += nanoTime() - t1;
    }

    public void doStepInternal() {
        if (tNotStep != 0) tNotStep += nanoTime();
        long t1 = nanoTime();
        super.doStepInternal();
//        System.out.println("step " + getStepCount());

        collisionUpdateCountdown--;
        if (collisionUpdateCountdown == 0) {
            double n = 1 + box.getLeafList().size();
            boolean forcedUpdates = neighborManager instanceof NeighborManager.NeighborEventSource;
            //System.out.println(intervalCollisionCount+" "+n);
            if (forcedUpdates || intervalCollisionCount / n < 20.0) {
                collisionUpdateInterval++;
                if (verbose) System.out.println("interval => " + collisionUpdateInterval);
                tMax = collisionUpdateInterval * timeStep;
                if (verbose) System.out.println(tMax / eventBinsFirstAtom.length + " time per bin");
            }
            intervalCollisionCount = 0;
            computeAllCollisions();
        }

        IAtomList atoms = box.getLeafList();
        while (true) {
            int c = nextCollider();
            nextColliderIndex = c;
            if (c < 0) {
                // we have no record of the next collision
                break;
            }
            double tcol = collisionTimes[c];
            if (tcol > tBase + timeStep) break;
            if (tcol < tBase - 1e-12) {
                throw new RuntimeException("nope " + c + " " + collisionPartners[c] + " " + tcol + " " + tBase);
            }
            int cPartner = collisionPartners[c];
            int oldState = collisionOldState[c];
            double[] virial = {0}, du = {0};
            long t1b = nanoTime();
            IAtomKinetic cAtom = (IAtomKinetic) atoms.get(c);
            if (cPartner >= 0) {
                IPotential2 pHard = collisionPotentials[c];
                IAtomKinetic partnerAtom = (IAtomKinetic) atoms.get(cPartner);
                Vector dr = collisionVector[c];
                Vector dv = space.makeVector();
                dv.Ev1Mv2(partnerAtom.getVelocity(), cAtom.getVelocity());

                int newState = pHard.bump(cAtom, partnerAtom, oldState, dr, dv, tcol - tBase, virial, du);
//                System.out.println("bump "+c+" "+cPartner+" at "+tcol+" state "+oldState+" => "+newState);
//                System.out.println(" dr "+dr+" "+Math.sqrt(dr.squared()));
                if (newState != oldState) {
                    setPairState(cAtom, partnerAtom, newState);
                }
                tBump += nanoTime() - t1b;

                double dE = du[0];
                currentPotentialEnergy += dE;
                currentKineticEnergy -= dE;
//                System.out.printf("%6d %3d %3d %7.3e %3.2f %2d %2d % f % f\n", collisionCount, Math.min(c, cPartner), Math.max(c, cPartner), tcol, Math.sqrt(dr.squared()), oldState, newState, dr.dot(dv), virial[0] / (dr.dot(dv)));
//                System.out.printf("%6d %3d %3d %3.1f % f % f\n", collisionCount, Math.min(c,cPartner), Math.max(c,cPartner), Math.sqrt(dr.squared()), dr.dot(dv), virial[0]/(dr.dot(dv)));

                long t1data = nanoTime();
                for (IntegratorHard.CollisionListener listener : this.collisionListeners) {
                    listener.pairCollision(cAtom, partnerAtom, dr, dv, virial[0], tcol - tBase);
                }
                tData += nanoTime() - t1data;
                updatePair(cAtom, partnerAtom, tcol - tBase);
            } else if (cPartner == -1) {
                IPotential1 pHard = fieldPotentials[cAtom.getType().getIndex()];
                Vector r = collisionVector[c];
                Vector deltaP = space.makeVector();
                fieldState[c] = pHard.bump(cAtom, fieldState[c], collisionVector[c], tcol - tBase, deltaP, du);
//                System.out.println("bump "+c+" at "+tcol+" state "+oldState+" => "+fieldState[c]);
                tBump += nanoTime() - t1b;

                double dE = du[0];
                currentPotentialEnergy += dE;
                currentKineticEnergy -= dE;
//                System.out.printf("%6d %3d %7.3e %2d %2d\n", collisionCount, c, tcol, oldState, fieldState[c]);
//            System.out.printf("%6d %3d %3d %3.1f % f % f\n", collisionCount, Math.min(c,cPartner), Math.max(c,cPartner), Math.sqrt(dr.squared()), dr.dot(dv), virial[0]/(dr.dot(dv)));

                long t1data = nanoTime();
                for (IntegratorHard.CollisionListener listener : this.collisionListeners) {
                    listener.fieldCollision(cAtom, r, deltaP, tcol - tBase);
                }
                tData += nanoTime() - t1data;
                updateAtom(cAtom, tcol - tBase);
            } else if (cPartner == -2) {
                // null potential collision
//                System.out.printf("%6d %3d %7.3e (null)\n", collisionCount, c, tcol);
                updateAtom(cAtom, tcol - tBase);
            }
            collisionCount++;
            intervalCollisionCount++;
        }
        tBase += timeStep;
//        System.out.println("advancing "+timeStep+" for "+stepCount);
        long t1a = nanoTime();
        advanceAcrossTimeStep(timeStep);

        if (isothermal) {
            doThermostatInternal();
        }

        if (writeTiming) {
            tAdvance += nanoTime() - t1a;
            tSteps += nanoTime() - t1;
            int s = box.getLeafList().size();
            if (s == 0) s = 1;
            int printInterval = 4000000 / s;
            if (stepCount >= lastTimingPrint + printInterval) {
                double scale = box.getLeafList().size() * printInterval;
                System.out.printf("time/step: %3.1f %3.1f  update: %3.1f  collect: %3.1f  up: %3.1f  down: %3.1f  ctime: %3.1f  add: %3.1f  next: %3.1f  delete: %3.1f  bump: %3.1f  data: %3.1f adv: %3.1f  nPerBin: %d  nPerBinMax: %d\n",
                        tSteps / scale, tNotStep / scale, tUpdate / scale, tCollect / scale, tUp / scale, tDown / scale, tCollision / scale, tAdd / scale, tNext / scale, tDelete / scale, tBump / scale, tData / scale, tAdvance / scale, (int) (nPerBinSum / (double) nComputeAll), nPerBinMax);
                tNotStep = tAdvance = tCollision = tUpdate = tCollect = tUp = tDown = tAdd = tSteps = tNext = tDelete = tBump = tData = 0;
                nPerBinMax = nPerBinSum = nComputeAll = 0;
                lastTimingPrint = stepCount;
            }
            tNotStep -= nanoTime();
        }

        // we may have advanced our bin into the next timestep
        lastBin = (int) (tBase / tMax * eventBinsFirstAtom.length);
    }

    /**
     * Updates collision times appropriately after randomizing momenta
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomenta() {
        super.randomizeMomenta();
        // super.randomizeMomenta alters the velocities, so we need to
        // recalculate collision times
        computeAllCollisions();
    }

    public void boxMoleculeAdded(BoxMoleculeEvent e) {
        // hook in here so that we can suppress computing a collision time now.
        // We'll compute the collision times later when the molecule's atoms get new
        // velocities.
        handlingEvent++;
        super.boxMoleculeAdded(e);
        handlingEvent--;
    }

    public void postRestore() {
        super.postRestore();
        computeAllCollisions();
    }

    protected void scaleMomenta() {
        super.scaleMomenta();
        // super.scaleMomenta alters the velocities, so we need to
        // recalculate collision times
        if (initialized) computeAllCollisions();
    }

    protected void randomizeTotalKE() {
        super.randomizeTotalKE();
        // super.randomizeMomenta alters the velocities, so we need to
        // recalculate collision times
        if (initialized) computeAllCollisions();
    }

    protected void randomizeMomentum(IAtomKinetic atom) {
        super.randomizeMomentum(atom);
        if (initialized && handlingEvent == 0) {
            updateAtom(atom, 0);
        }
    }

    protected void setPairState(IAtom iAtom, IAtom jAtom, int newState) {
        int i = iAtom.getLeafIndex();
        IMolecule parentMolecule = iAtom.getParentGroup();
        if (bondingInfo != null && jAtom.getParentGroup() == parentMolecule) {
            boolean[] found = {false};
            bondingInfo.bondedPairPartners[iAtom.getParentGroup().getType().getIndex()].forEach((potential, partners) -> {
                for (int[] pairIdx : partners[iAtom.getIndex()]) {
                    IAtom atom1 = parentMolecule.getChildList().get(pairIdx[0]);
                    IAtom atom2 = parentMolecule.getChildList().get(pairIdx[1]);
                    IAtom jjAtom = atom1 == iAtom ? atom2 : atom1;
                    if (jjAtom.getLeafIndex() < i) continue;
                    if (jjAtom != jAtom) {
                        continue;
                    }
                    bondState[i].put(jAtom.getLeafIndex(), newState);
                    found[0] = true;
                    return;
                }
            });
            if (found[0]) return;
        }
        neighborManager.setPairState(i, jAtom.getLeafIndex(), newState);
    }

    /**
     * Advances all atom coordinates by tStep, without any intervening collisions.
     * Uses free-flight kinematics.
     */
    protected void advanceAcrossTimeStep(double tStep) {
        IAtomList leafList = box.getLeafList();
        for (IAtom iAtom : leafList) {
            iAtom.getPosition().PEa1Tv1(tStep, ((IAtomKinetic) iAtom).getVelocity());
        }
    }

    public IAtom[] getNextCollidingPair() {
        if(nextColliderIndex < 0) {
            nextCollidingPair[0] = null;
            nextCollidingPair[1] = null;
        } else {
            nextCollidingPair[0] = box.getLeafList().get(nextColliderIndex);
            int partner = collisionPartners[nextColliderIndex];
            nextCollidingPair[1] = (partner >= 0) ? box.getLeafList().get(partner) : null;
        }
        return nextCollidingPair;
    }

    @Override
    public void neighborListNeighborsUpdated() {
        computeAllCollisions();
    }

    public interface CollisionListener {
        void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial, double tCollision);

        default void fieldCollision(IAtomKinetic atom, Vector r, Vector deltaP, double tCollision) {
        }
    }

    // called by neighbor consumer and also bonded callback
    private interface MyConsumer {
        void accept(IAtom jAtom, Vector rij, int state, IPotential2 pij);
    }
}

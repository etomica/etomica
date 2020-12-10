/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.integrator;

import etomica.atom.*;
import etomica.box.Box;
import etomica.nbr.list.INeighborListListener;
import etomica.nbr.list.NeighborListManagerFasterer;
import etomica.potential.IPotentialHard;
import etomica.potential.Potential2Soft;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialComputePair;
import etomica.space.Vector;
import etomica.util.collections.Int2IntHash;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IntegratorHardFasterer extends IntegratorMDFasterer implements INeighborListListener {

    protected final NeighborManager neighborManager;
    protected double[] collisionTimes;
    protected int[] collisionPartners;
    protected IPotentialHard[] collisionPotentials;
    protected Vector[] collisionVector;
    protected int[] eventBinsFirstAtom, eventBinsNextAtom, eventBinsPrevAtom;
    protected final IPotentialHard[][] pairPotentials;
    protected double tBase, tMax;
    protected int collisionUpdateInterval, collisionUpdateCountdown;
    protected long collisionCount;
    protected final NeighborIterator neighborIterator;
    protected final AtomLeafAgentManager<Int2IntHash> stateManager;
    protected final List<CollisionListener> collisionListeners;
    private int lastBin = 0;

    private long tDelete, tAdd, tNext, tSteps, tUpdate, tUp, tDown, tBump, tData, tCollect, tCollision, tAdvance, tNotStep;
    private final boolean writeTiming = false, verbose = false;

    /**
     * Constructs integrator with a default for non-isothermal sampling.
     *
     * @param potentialCompute PotentialMaster instance used to compute the energy and forces
     * @param random           random number generator used for initial velocities and some thermostats
     * @param timeStep         time step for integration
     * @param temperature      used by thermostat and/or to initialize velocities
     * @param box
     */
    public IntegratorHardFasterer(PotentialComputePair potentialCompute, NeighborManager neighborManager, IRandom random, double timeStep, double temperature, Box box) {
        super(potentialCompute, random, timeStep, temperature, box);
        this.neighborManager = neighborManager;
        eventBinsFirstAtom = new int[10];
        collisionVector = new Vector[0];
        resizeArrays();
        int n = potentialCompute.getPairPotentials().length;
        pairPotentials = new IPotentialHard[n][n];
        // we want to store enough collisions that we don't have to recompute them
        // before the neighbor update
        collisionUpdateInterval = 5;
        tMax = collisionUpdateInterval * timeStep;
        neighborIterator = neighborManager.makeNeighborIterator();
        stateManager = new AtomLeafAgentManager<>(iAtom -> new Int2IntHash(10), box);
        if (neighborManager instanceof NeighborListManagerFasterer) {
            ((NeighborListManagerFasterer) neighborManager).addListener(this);
        }
        collisionListeners = new ArrayList<>(0);
    }

    public void setCollisionUpdateInterval(int newInterval) {
        collisionUpdateInterval = newInterval;
        tMax = collisionUpdateInterval * timeStep;
        if (initialized) reset();
    }

    public void addCollisionListener(CollisionListener newListener) {
        collisionListeners.add(newListener);
    }

    public void setTimeStep(double dt) {
        super.setTimeStep(dt);
        tMax = collisionUpdateInterval * timeStep;
        if (initialized) reset();
    }

    protected void resizeArrays() {
        int nAtoms = box.getLeafList().size();
        if (collisionVector.length < nAtoms) {
            collisionTimes = new double[nAtoms];
            collisionPartners = new int[nAtoms];
            eventBinsNextAtom = new int[nAtoms];
            eventBinsPrevAtom = new int[nAtoms];
            collisionPotentials = new IPotentialHard[nAtoms];
            int old = collisionVector.length;
            collisionVector = new Vector[nAtoms];
            for (int i = old; i < nAtoms; i++) {
                collisionVector[i] = space.makeVector();
            }
        }
    }

    public void reset() {
        super.reset();
        resizeArrays();
        Potential2Soft[][] softPotentials = ((PotentialComputePair) potentialCompute).getPairPotentials();
        int n = softPotentials.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                pairPotentials[i][j] = (IPotentialHard) softPotentials[i][j];
            }
        }
        neighborManager.init();

        collisionUpdateCountdown = 1;
        IAtomList atoms = box.getLeafList();
        for (IAtom atom : atoms) {
            IAtomKinetic iAtom = (IAtomKinetic) atom;
            int iType = iAtom.getType().getIndex();
            Int2IntHash hash = stateManager.getAgent(iAtom);
            hash.clear();

            neighborIterator.iterDownNeighbors(iAtom.getLeafIndex(), (jAtom, rij) -> {
                int j = jAtom.getLeafIndex();
                IAtomKinetic jjAtom = (IAtomKinetic) jAtom;
                int state = pairPotentials[iType][jAtom.getType().getIndex()].getState(iAtom, jjAtom, rij);
                if (state > -1) hash.put(j, state);
            });
        }
        computeAllCollisions();
    }

    protected long nPerBinSum = 0, nComputeAll = 0, nPerBinMax = 0;

    protected void computeAllCollisions() {
        tBase = 0;
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
            eventBinsFirstAtom = new int[eventBinsFirstAtom.length * 2];
            if (verbose) System.out.println("bins increased to " + eventBinsFirstAtom.length);
            Arrays.fill(eventBinsFirstAtom, -1);
            IAtomList atoms = box.getLeafList();
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
        long t1 = System.nanoTime();
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
                tNext += System.nanoTime() - t1;
                return first;
            }
        }
        lastBin = eventBinsFirstAtom.length;
        tNext += System.nanoTime() - t1;
        return -1;
    }

    protected void collectDownColliders(IAtom atom, AtomArrayList downColliders) {
        int i = atom.getLeafIndex();
        neighborIterator.iterDownNeighbors(atom.getLeafIndex(), (jAtom, rij) -> {
            if (collisionPartners[jAtom.getLeafIndex()] == i) {
                downColliders.add(jAtom);
            }
        });
    }

    protected void removeCollision(int i) {
        long t1 = System.nanoTime();
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
        tDelete += System.nanoTime() - t1;
    }

    protected void collisionTimeDown(IAtomKinetic iAtom, double falseTime) {
        long t1cd = System.nanoTime();
        Vector iVelocity = iAtom.getVelocity();
        int iType = iAtom.getType().getIndex();
        int i = iAtom.getLeafIndex();
        neighborIterator.iterDownNeighbors(iAtom.getLeafIndex(), (jAtom, rij) -> {
            int j = jAtom.getLeafIndex();
            IAtomKinetic jjAtom = (IAtomKinetic) jAtom;
            Vector dv = space.makeVector();
            dv.Ev1Mv2(jjAtom.getVelocity(), iVelocity);
            rij.PEa1Tv1(falseTime, dv);
            long t1 = System.nanoTime();
            Int2IntHash jHash = stateManager.getAgent(jAtom);
            double time = pairPotentials[iType][jAtom.getType().getIndex()].collisionTime(iAtom, jjAtom, rij, dv, jHash.get(i));
            tCollision += System.nanoTime() - t1;
            if (time < 0) {
                System.out.println("computed col down " + i + " with " + j + " at " + time);
                Vector dr = space.makeVector();
                dr.Ev1Mv2(jAtom.getPosition(), iAtom.getPosition());
                System.out.println("dr " + dr + " " + dr.squared());
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
                collisionPotentials[j] = pairPotentials[iType][jAtom.getType().getIndex()];
                rij.PEa1Tv1(time, dv);
                collisionVector[j].Ea1Tv1(-1, rij);

                long t1a = System.nanoTime();
                int b = (int) (collisionTimes[j] / tMax * eventBinsFirstAtom.length);

                int oldFirst = eventBinsFirstAtom[b];
                eventBinsNextAtom[j] = oldFirst;
                eventBinsPrevAtom[j] = -1 - b;
                if (oldFirst >= 0) eventBinsPrevAtom[oldFirst] = j;
                eventBinsFirstAtom[b] = j;
                tAdd += System.nanoTime() - t1a;
            }
        });

        tDown += System.nanoTime() - t1cd;
    }

    protected void collisionTimeUp(IAtomKinetic iAtom, double falseTime) {
        long t1cu = System.nanoTime();
        int i = iAtom.getLeafIndex();
        if (collisionTimes[i] < tMax) {
            // need to remove i from event bins
            removeCollision(i);
        }
        collisionTimes[i] = Double.POSITIVE_INFINITY;
        Vector iVelocity = iAtom.getVelocity();
        int iType = iAtom.getType().getIndex();
        double[] tCheck = new double[]{tMax - tBase - falseTime};
        Int2IntHash iHash = stateManager.getAgent(iAtom);
        neighborIterator.iterUpNeighbors(iAtom.getLeafIndex(), (jAtom, rij) -> {

            IAtomKinetic jjAtom = (IAtomKinetic) jAtom;
            Vector dv = space.makeVector();
            dv.Ev1Mv2(jjAtom.getVelocity(), iVelocity);
            rij.PEa1Tv1(falseTime, dv);
            long t1ct = System.nanoTime();
            int j = jAtom.getLeafIndex();
            double time = pairPotentials[iType][jAtom.getType().getIndex()].collisionTime(iAtom, jjAtom, rij, dv, iHash.get(j));
            tCollision += System.nanoTime() - t1ct;
            if (time < 0) {
                throw new RuntimeException("Negative up collision time between " + i + " and " + jAtom.getLeafIndex() + " at " + time);
            }
            if (time < tCheck[0]) {
//                long t1 = System.nanoTime();
                collisionTimes[i] = time;
                tCheck[0] = time;
                collisionPartners[i] = j;
                collisionPotentials[i] = pairPotentials[iType][jAtom.getType().getIndex()];
                rij.PEa1Tv1(time, dv);
                collisionVector[i].E(rij);
//                tAdd += System.nanoTime() - t1;
            }
        });

        collisionTimes[i] += tBase + falseTime;
        long t1 = System.nanoTime();
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
        tAdd += System.nanoTime() - t1;
        tUp += System.nanoTime() - t1cu;
    }

    public void updatePair(IAtomKinetic atom1, IAtomKinetic atom2, double falseTime) {
        long t1 = System.nanoTime();
        AtomArrayList downColliders = new AtomArrayList(0);
        collectDownColliders(atom1, downColliders);
        collectDownColliders(atom2, downColliders);
        tCollect += System.nanoTime() - t1;
        for (IAtom iAtom : downColliders) {
            collisionTimeUp((IAtomKinetic) iAtom, falseTime);
        }

        // don't need to do up for atom1 because it was included above!
        collisionTimeDown(atom1, falseTime);

        collisionTimeUp(atom2, falseTime);
        collisionTimeDown(atom2, falseTime);
        tUpdate += System.nanoTime() - t1;
    }

    public void doStepInternal() {
        if (tNotStep != 0) tNotStep += System.nanoTime();
        long t1 = System.nanoTime();
        super.doStepInternal();

        collisionUpdateCountdown--;
        if (collisionUpdateCountdown == 0) {
            collisionUpdateInterval++;
            if (verbose) System.out.println("interval => " + collisionUpdateInterval);
            tMax = collisionUpdateInterval * timeStep;
            computeAllCollisions();
        }

        IAtomList atoms = box.getLeafList();
        while (true) {
            int c = nextCollider();
            if (c < 0) {
                // we have no record of the next collision
                break;
            }
            double tcol = collisionTimes[c];
            if (tcol > tBase + timeStep) break;
            int cPartner = collisionPartners[c];
            IPotentialHard pHard = collisionPotentials[c];
            IAtomKinetic cAtom = (IAtomKinetic) atoms.get(c);
            IAtomKinetic partnerAtom = (IAtomKinetic) atoms.get(cPartner);
            Vector dr = collisionVector[c];
            Vector dv = space.makeVector();
            dv.Ev1Mv2(partnerAtom.getVelocity(), cAtom.getVelocity());

            if (tcol < tBase - 1e-12) {
                throw new RuntimeException("nope");
            }
            Int2IntHash agent = stateManager.getAgent(cAtom);
            int oldState = agent.get(cPartner);
            double[] virial = {0}, du = {0};
            long t1b = System.nanoTime();
            int newState = pHard.bump(cAtom, partnerAtom, oldState, dr, dv, tcol - tBase, virial, du);
            if (newState < 0) {
                if (oldState >= 0) {
                    agent.remove(cPartner);
                }
            } else if (newState != oldState) {
                agent.put(cPartner, newState);
            }
            tBump += System.nanoTime() - t1b;
            long t1data = System.nanoTime();
            double dE = du[0];
            currentPotentialEnergy += dE;
            currentKineticEnergy -= dE;

            for (IntegratorHardFasterer.CollisionListener listener : this.collisionListeners) {
                listener.collisionAction(cAtom, partnerAtom, dr, dv, virial[0]);
            }
            collisionCount++;
            tData += System.nanoTime() - t1data;
            updatePair(cAtom, partnerAtom, tcol - tBase);
        }
        tBase += timeStep;
//        System.out.println("advancing "+timeStep);
        long t1a = System.nanoTime();
        advanceAcrossTimeStep(timeStep);
        tAdvance += System.nanoTime() - t1a;
        tSteps += System.nanoTime() - t1;
        int printInterval = 4000000 / box.getLeafList().size();
        if (writeTiming && stepCount % printInterval == 0) {
            double scale = box.getLeafList().size() * printInterval;
            System.out.printf("time/step: %3.1f %3.1f  update: %3.1f  collect: %3.1f  up: %3.1f  down: %3.1f  ctime: %3.1f  add: %3.1f  next: %3.1f  delete: %3.1f  bump: %3.1f  data: %3.1f adv: %3.1f  nPerBin: %d  nPerBinMax: %d\n",
                    tSteps / scale, tNotStep / scale, tUpdate / scale, tCollect / scale, tUp / scale, tDown / scale, tCollision / scale, tAdd / scale, tNext / scale, tDelete / scale, tBump / scale, tData / scale, tAdvance / scale, (int) (nPerBinSum / (double) nComputeAll), nPerBinMax);
            tNotStep = tAdvance = tCollision = tUpdate = tCollect = tUp = tDown = tAdd = tSteps = tNext = tDelete = tBump = tData = 0;
            nPerBinMax = nPerBinSum = nComputeAll = 0;
        }
        tNotStep -= System.nanoTime();
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

    @Override
    public void neighborListNeighborsUpdated() {
        computeAllCollisions();
    }

    public interface CollisionListener {
        void collisionAction(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial);
    }
}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.IListener;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;

/**
 * Looks for atoms without a full set of first nearest neighbors and attempts
 * insertions adjacent to one of those atoms.  Hoops are jumped through to
 * ensure detailed balance.
 * 
 * Neighbor lists are used to find insertion points and candidates for
 * deletion.  Internal lists are constructed (relatively expensive) and then
 * reconstructed when a move is accepted or when some other part of the
 * simulation moves atoms (this part is somewhat hardcoded; enjoy).
 *
 * @author Andrew Schultz
 */
public class MCMoveInsertDeleteLatticeVacancy extends MCMoveInsertDeleteBiased implements IListener<MCMoveEvent> {

    protected final NeighborIterator neighborIterator;
    protected final Vector dest;
    protected Integrator integrator;
    protected long lastStepCount;
    protected boolean dirty;
    protected double maxDistance, maxInsertDistance;
    protected double nbrDistance;
    protected IntArrayList insertCandidates, deleteCandidates;
    protected int[] numNeighbors, numNeighborCandidatesOnDelete, deleteCandidateTimes, numDeleteCandidateNbrs;
    protected int totalDeleteCandidateTimes;
    protected int numNewDeleteCandidates;
    protected int forced = 0;
    protected double oldLnA, oldB;
    protected final Vector oldPosition;
    protected Vector[] nbrVectors;
    protected final Space space;
    protected double oldBoxSize;

    public MCMoveInsertDeleteLatticeVacancy(PotentialCompute potentialMaster, NeighborManager neighborManager,
                                            IRandom random, Space _space, Integrator integrator, double maxDistance, int maxN, int maxVacancy) {
        super(potentialMaster, random, _space, maxN-maxVacancy, maxN);
        this.space = _space;
        this.neighborIterator = neighborManager.makeNeighborIterator();
        dest = _space.makeVector();
        this.integrator = integrator;
        this.maxDistance = maxDistance;
        maxInsertDistance = 0.05;
        insertCandidates = new IntArrayList();
        deleteCandidates = new IntArrayList();
        numNeighbors = new int[0];
        numNeighborCandidatesOnDelete = new int[0];
        lastStepCount = -1;
        dirty = true;
        lnbias = new double[0];
        oldPosition = _space.makeVector();
    }

    public void setBox(final Box box) {
        super.setBox(box);
        oldBoxSize = box.getBoundary().getBoxSize().getX(0);
    }

    public void setFluidNbrDistance(double newNbrDistance) {
        nbrDistance = newNbrDistance;
        double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
        if (maxInsertNbrDistance > maxDistance) {
            throw new RuntimeException("nbrDistance must be greater than maxInsert distance");
        }
        nbrVectors = null;
    }

    public void makeFccVectors(double newNbrDistance) {
        nbrDistance = newNbrDistance;
        double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
        if (maxInsertNbrDistance > maxDistance) {
            throw new RuntimeException("nbrDistance must be greater than maxInsert distance");
        }
        nbrVectors = new Vector[12];
        double s = nbrDistance/Math.sqrt(2);
        for (int i=0; i<12; i++) {
            nbrVectors[i] = space.makeVector();
            boolean even = i%2 == 0;

            if (i < 4) {
                nbrVectors[i].setX(0, i<2 ? -s : s);
                nbrVectors[i].setX(1, even ? -s : s);
            }
            else if (i < 8) {
                nbrVectors[i].setX(1, i<6 ? -s : s);
                nbrVectors[i].setX(2, even ? -s : s);
            }
            else {
                nbrVectors[i].setX(0, i<10 ? -s : s);
                nbrVectors[i].setX(2, even ? -s : s);
            }
        }
    }

    public void makeHcpVectors(double newNbrDistance) {
        nbrDistance = newNbrDistance;
        double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
        if (maxInsertNbrDistance > maxDistance) {
            throw new RuntimeException("nbrDistance must be greater than maxInsert distance");
        }
        nbrVectors = new Vector[12];
        double s = nbrDistance;

        for (int i=0; i<12; i++) {
            nbrVectors[i] = space.makeVector();

            if (i < 6) {
                nbrVectors[i].setX(0, s*Math.cos(i*Math.PI/3));
                nbrVectors[i].setX(1, s*Math.sin(i*Math.PI/3));
            }
            else {
                nbrVectors[i].setX(2, 0.5*s*(i<9?(-1):(+1))*Math.sqrt(8.0/3.0));
                // 0.5*s, 0.5*sqrt(3)/2*s
                // 0.25*s^2 + 0.25*3/4 = 0.25*s^2*7/4
                nbrVectors[i].setX(0, s*Math.cos((2*i + 0.5)*Math.PI/3)/Math.sqrt(3));
                nbrVectors[i].setX(1, s*Math.sin((2*i + 0.5)*Math.PI/3)/Math.sqrt(3));
            }
        }
    }

    public void setMaxDistance(double maxDistance) {
        this.maxDistance = maxDistance;
    }

    public double getMaxDistance() {
        return maxDistance;
    }

    public void setMaxInsertDistance(double maxInsertDistance) {
        this.maxInsertDistance = maxInsertDistance;
    }

    public void reset() {
        dirty = true;
    }

    public void setLnBias(int n, double nBias) {
        forced = 0;
        super.setLnBias(n, nBias);
    }

    public boolean doTrial() {
        if (integrator.getStepCount() != lastStepCount) {
            forced = 0;
        }
        int numAtoms = box.getLeafList().size();
        if (forced==1) {
            insert = !insert;
            System.out.println("forcing "+(insert ? "insertion" : "deletion"));
            forced=2;
        }
        else {
            insert = (random.nextInt(2) == 0);
        }
        if (insert && numAtoms == maxN) return false;
        if (!insert && numAtoms == minN) return false;
        numNewDeleteCandidates = 0;
        if (dirty || lastStepCount != integrator.getStepCount()) findCandidates();
        if (insert) {
            if(!reservoir.isEmpty()) testMolecule = reservoir.remove(reservoir.size()-1);
            else testMolecule = species.makeMolecule();
            IAtom testAtom = testMolecule.getChildList().get(0);

            int nInsertCandidates = insertCandidates.size();
            if (nInsertCandidates == 0) {
                reservoir.add(testMolecule);
                // test molecule is no longer in the simulation and should not be 
                // returned by affectedAtoms
                testMolecule = null;
                return false;
            }
            IAtom partner = box.getLeafList().get(insertCandidates.getInt(random.nextInt(nInsertCandidates)));
//            System.out.println("inserting next to "+partner);
            uOld = 0;
            if (nbrVectors != null) {
                int iLatVec = random.nextInt(nbrVectors.length);
                dest.setRandomInSphere(random);
                dest.TE(maxInsertDistance);
                dest.PE(nbrVectors[iLatVec]);
            }
            else {
                //XXX still wrong
                dest.setRandomInSphere(random);
                dest.TE(maxInsertDistance);
                Vector dr = space.makeVector();
                dr.E(dest);
                dest.setRandomSphere(random);
                dest.TE(nbrDistance);
                dest.PE(dr);
            }
            testAtom.getPosition().E(partner.getPosition());
            testAtom.getPosition().PE(dest);
            testAtom.getPosition().PE(box.getBoundary().centralImage(testAtom.getPosition()));
            if (forced==2) {
                testAtom.getPosition().E(oldPosition);
            }
            box.addMolecule(testMolecule);
            uNew = potentialMaster.computeOne(testAtom);

            // inserting testAtom might change some deleteCandidates
            // some existing deleteCandidates with nbrs=12 might have 13 now
            // we also need to see how many times testAtom shows up as a neighbor
            // of a deleteCandidate

            int[] nTestNbrs = {0}, nTestNbrsDeletion = {0};
            double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
            double minInsertNbrDistance = nbrDistance - maxInsertDistance;
            neighborIterator.iterAllNeighbors(testAtom.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij, int n) {
                    int jj = jAtom.getLeafIndex();
                    double r2 = rij.squared();
                    if (r2 < maxDistance*maxDistance) {
                        nTestNbrs[0]++;
                        boolean deleteCand = false;
                        if (r2 < maxInsertNbrDistance*maxInsertNbrDistance && r2 > minInsertNbrDistance*minInsertNbrDistance) {
                            for (int k=0; k<nbrVectors.length; k++) {
                                double s2 = rij.Mv1Squared(nbrVectors[k]);
                                if (s2 < maxInsertDistance*maxInsertDistance) {
                                    deleteCand = true;
                                    break;
                                }
                            }
                            if (deleteCand) {
                                nTestNbrsDeletion[0]++;
                            }
                        }
                        if (numNeighbors[jj] == 12) {
                            // jj now has 13 neighbors
                            // by inserting testAtom, jj's old neighbors can no longer be deleted
                            numNewDeleteCandidates-=numDeleteCandidateNbrs[jj];
                        }
                        else if (numNeighbors[jj] < 12 && deleteCand) {
                            // testAtom would be a deleteCandidate due to jj
                            numNewDeleteCandidates++;
                        }
                    }
                }
            });

            // it could be that testAtom itself has 12 or fewer neighbors
            // if so, those neighbors would now be delete candidates
            if (nTestNbrs[0]<=12) {
                numNewDeleteCandidates += nTestNbrsDeletion[0];
                
            }
        } else {//delete
            if(box.getNMolecules(species) == 0) {
                testMolecule = null;
                return false;
            }
            int nDeleteCandidates = deleteCandidates.size();
            if (nDeleteCandidates == 0) return false;
            int irand = random.nextInt(totalDeleteCandidateTimes);
            int icount = 0;
            int ip = -1;
            for (int ii=0; ii<deleteCandidates.size(); ii++) {
                int i = deleteCandidates.getInt(ii);
                icount += deleteCandidateTimes[i];
                if (icount > irand) {
                    ip = i;
                    break;
                }
            }
            if (forced==2) {
                ip = numAtoms-1;
            }

            IAtom testAtom = box.getLeafList().get(ip);
            testMolecule = testAtom.getParentGroup();
            //delete molecule only upon accepting trial
            uOld = potentialMaster.computeOneOld(testAtom);
            
            // by deleting testAtom, we might turn other atoms into "insertCandidates"
            // numNeighborCandidates is how many of our neighbors have 12 nbrs.
            uNew = 0;
            
        }
        return true;
    }

    protected void findCandidates() {
        double newBoxSize = box.getBoundary().getBoxSize().getX(0);
        if (Math.abs(newBoxSize-oldBoxSize)/oldBoxSize > 1e-14) {
            for (int i=0; i<12; i++) {
                nbrVectors[i].TE(newBoxSize/oldBoxSize);
            }
            nbrDistance *= newBoxSize/oldBoxSize;
            oldBoxSize = newBoxSize;
        }
        
        int numAtoms = box.getLeafList().size();
        if (numNeighbors.length < numAtoms) {
            numNeighbors = new int[numAtoms];
            numNeighborCandidatesOnDelete = new int[numAtoms];
            numDeleteCandidateNbrs = new int[numAtoms];
            deleteCandidateTimes = new int[numAtoms];
        }
        insertCandidates.clear();
        deleteCandidates.clear();
        totalDeleteCandidateTimes = 0;
        for (int i=0; i<numAtoms; i++) {
            numNeighbors[i] = numNeighborCandidatesOnDelete[i] = numDeleteCandidateNbrs[i] = deleteCandidateTimes[i] = 0;
        }
        for (int i=0; i<numAtoms; i++) {
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij, int n) {
                    double r2 = rij.squared();
                    if (r2 < maxDistance*maxDistance) {
                        numNeighbors[finalI]++;
                        numNeighbors[jAtom.getLeafIndex()]++;
                    }
                }
            });
        }
        double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
        double minInsertNbrDistance = nbrDistance - maxInsertDistance;
        for (int i=0; i<numAtoms; i++) {
            if (numNeighbors[i] < 13) {
                // the neighbors of i may be candidates for deletion.  after deleting
                // one of its neighbors, i would have <12 neighbors
                int finalI = i;
                neighborIterator.iterAllNeighbors(i, new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom jAtom, Vector rij, int n) {
                        int jj = jAtom.getLeafIndex();
                        double r2 = rij.squared();
                        if (numNeighbors[finalI] == 12 && r2 < maxDistance*maxDistance) {
                            // if we delete jj then i becomes an "insertCandidate"
                            numNeighborCandidatesOnDelete[jj]++;
                        }
                        if (r2 < maxInsertNbrDistance*maxInsertNbrDistance && r2 > minInsertNbrDistance*minInsertNbrDistance) {
                            boolean success = false;
                            for (int k=0; k<nbrVectors.length; k++) {
                                double s2 = rij.Mv1Squared(nbrVectors[k]);
                                if (s2 < maxInsertDistance*maxInsertDistance) {
                                    success = true;
                                    break;
                                }
                            }
                            if (success) {
                                // we need to know how many times jj shows up as a delete candidate
                                deleteCandidateTimes[jj]++;
                                numDeleteCandidateNbrs[finalI]++;
                                totalDeleteCandidateTimes++;
                                if (deleteCandidateTimes[jj] > 1) return;
                                deleteCandidates.add(jj);
                            }
                        }
                    }
                });
                if (numNeighbors[i] < 12) {
                    // we will attempt to insert next to i
                    insertCandidates.add(i);
                    // deleting i will remove an insert candidate
                    numNeighborCandidatesOnDelete[i]--;
                }
            }
        }
        dirty = false;
        lastStepCount = integrator.getStepCount();
    }

    public double getChi(double temperature) {
        double lna = getLnBiasDiff();

        double shellV = nbrVectors.length*4.0/3.0*Math.PI*maxInsertDistance*maxInsertDistance*maxInsertDistance;
        double c = 0;
        if (insert) {
            c = insertCandidates.size()*shellV/(totalDeleteCandidateTimes + numNewDeleteCandidates);
        }
        else {
            // our candidate for deletion was listed deleteCandidateTimes times
            int newInsertCandidates = insertCandidates.size() + numNeighborCandidatesOnDelete[testMolecule.getIndex()];
            c = totalDeleteCandidateTimes/(shellV*newInsertCandidates);
        }
        if (forced > 0) {
            double x = lna + Math.log(c);
            System.out.println("**** lnA = "+x+"    "+lna+"    "+Math.log(c));
            System.out.println("*old lnA = "+oldLnA);
            if (Math.abs(x+oldLnA) > 1e-6) {
                throw new RuntimeException("oops");
            }
        }
        oldLnA = lna+Math.log(c);
        if (false) {
            System.out.println(insert+" lnbias "+lna+" log(c) "+Math.log(c)+" log(a) "+oldLnA);
        }
        double A = Math.exp(lna) * c;

        double b = uOld - uNew;
        if (false) System.out.println(insert+" b = "+b);
        if (forced==2) {
            System.out.println("**** b = "+b);
            System.out.println("*old b = "+oldB);
            if (Math.abs((b+oldB)/(1+b-oldB)) > 1e-6) {
                throw new RuntimeException("oops");
            }
            forced = -1;
            return b;
        }
        else if (forced == -1) {
            forced = 0;
            return -Double.POSITIVE_INFINITY;
        }
        oldB = b;
        return A * Math.exp(b / temperature);
    }

    public void myAcceptNotify() {
        if (!insert) {
            oldPosition.E(testMolecule.getChildList().get(0).getPosition());
        }
        super.myAcceptNotify();
        dirty = true;
    }

    public void actionPerformed(MCMoveEvent event) {
        // any other move is accepted and we are dirty
        if (event instanceof MCMoveTrialCompletedEvent && event.getMCMove() != this && ((MCMoveTrialCompletedEvent)event).isAccepted()) {
            dirty = true;
        }
    }
}

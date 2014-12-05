/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import java.util.ArrayList;
import java.util.List;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IIntegrator;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.iterator.AtomIterator;
import etomica.nbr.list.NeighborListManager;
import etomica.space.ISpace;

/**
 * This move identifies vacancies as in the superclass, but then performs moves
 * which move an adjacent atom into a vacancy (# of atoms remains fixed).
 *  
 * @author Andrew Schultz
 */
public class MCMoveInsertDeleteLatticeVacancyJump extends
        MCMoveInsertDeleteLatticeVacancy {

    protected IAtom chosenAtom, testAtom;
    protected final AtomPair pair;
    protected final IPotentialAtomic p2;
    protected List<Integer>[] deleteCandidateNbrsOnDelete;
    
    public MCMoveInsertDeleteLatticeVacancyJump(
            IPotentialMaster potentialMaster, IRandom random, ISpace _space,
            IIntegrator integrator, double nbrDistance, IPotentialAtomic p2) {
        super(potentialMaster, random, _space, integrator, nbrDistance, 10, 0);
        pair = new AtomPair();
        this.p2 = p2;
    }
    
    public boolean doTrial() {
        numNewDeleteCandidates = 0;
        if (dirty || lastStepCount < integrator.getStepCount()) findCandidates();
        
        //delete
        int nDeleteCandidates = deleteCandidates.size();
        if (nDeleteCandidates == 0) return false;
        int irand = random.nextInt(totalDeleteCandidateTimes);
        int icount = 0;
        int ip = -1;
        for (int i : deleteCandidates) {
            icount += deleteCandidateTimes[i];
            if (icount > irand) {
                ip = i;
                break;
            }
        }
//        System.out.println("deleting "+ip);
        chosenAtom = box.getLeafList().getAtom(ip);
        //delete molecule only upon accepting trial
        energyMeter.setTarget(chosenAtom);
        uOld = energyMeter.getDataAsScalar();
        
        numNewDeleteCandidates = -1;
        // by deleting ip, each of these atoms' neighbors became eligible for deletion
        for (int jj : deleteCandidateNbrsOnDelete[ip]) {
            numNewDeleteCandidates += numDeleteCandidateNbrs[jj];
        }
        
        // by deleting testAtom, we might turn other atoms into "insertCandidates"
        // numNeighborCandidates is how many of our neighbors have 12 nbrs.

        // insert
        int nInsertCandidates = insertCandidates.size()+numNeighborCandidatesOnDelete[ip];
        if (nInsertCandidates == 0) return false;
        int ipc = random.nextInt(nInsertCandidates);
        if (ipc >= insertCandidates.size()) {
            // we've chosen to insert around an atom that was an insertion candidate
            // only because of the deleted atom (meaning that we insert near the deleted
            // atom).  We'll exclude this below, so just bail now.
            return false;
        }
        int ic = insertCandidates.get(ipc);
        IAtom partner = box.getLeafList().getAtom(ic);
//            System.out.println("inserting next to "+partner);

        int iLatVec = random.nextInt(nbrVectors.length);
        dest.setRandomInSphere(random);
        dest.TE(maxInsertDistance);
        dest.PE(nbrVectors[iLatVec]);
        if (testAtom == null) {
            testMolecule = species.makeMolecule();
            testAtom = testMolecule.getChildList().getAtom(0);
        }
        else {
            testMolecule = testAtom.getParentGroup();
        }
        
        testAtom.getPosition().E(partner.getPosition());
        testAtom.getPosition().PE(dest);
        if (testAtom.getPosition().Mv1Squared(chosenAtom.getPosition()) < maxDistance*maxDistance) {
            // insertion site and deleted atom are close enough that we would need to do
            // extra bookkeeping.  we could.  but it's easier to just bail and wait for
            // an easier jump to come along.
            return false;
        }

        box.addMolecule(testMolecule);
        
        pair.atom0 = chosenAtom;
        pair.atom1 = testAtom;
        uNew = -p2.energy(pair);
        if (uNew < -1e5) {
            // we inserted very close to chosen atom -- in the vacancy left by deleting it
            // we could fix this, but the whole move is pointless
            box.removeMolecule(testMolecule);
            return false;
        }
        // we inserted away from chosen atom
        energyMeter.setTarget(testAtom);
        uNew += energyMeter.getDataAsScalar();

        // inserting testAtom might change some deleteCandidates
        // some existing deleteCandidates with nbrs=12 might have 13 now
        // we also need to see how many times testAtom shows up as a neighbor
        // of a deleteCandidate

        IVector pi = testAtom.getPosition();
        IAtomList nbrs = potentialMaster.getNeighborManager(box).getUpList(testAtom)[0];
        int nTestNbrs = 0, nTestNbrsDeletion = 0;
        double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
        double minInsertNbrDistance = nbrDistance - maxInsertDistance;
        for (int j=0; j<nbrs.getAtomCount(); j++) {
            IAtom jAtom = nbrs.getAtom(j);
            int jj = jAtom.getLeafIndex();
            dr.Ev1Mv2(pi, jAtom.getPosition());
            box.getBoundary().nearestImage(dr);
            double r2 = dr.squared();
            if (r2 < maxDistance*maxDistance) {
                nTestNbrs++;
                boolean deleteCand = false;
                if (r2 < maxInsertNbrDistance*maxInsertNbrDistance && r2 > minInsertNbrDistance*minInsertNbrDistance) {
                    for (int k=0; k<nbrVectors.length; k++) {
                        double s2 = dr.Mv1Squared(nbrVectors[k]);
                        if (s2 < maxInsertDistance*maxInsertDistance) {
                            deleteCand = true;
                            break;
                        }
                    }
                    if (deleteCand) {
                        nTestNbrsDeletion++;
                    }
                }
                if (numNeighbors[jj] == 12) {
                    // jj now has 13 neighbors
                    // by inserting testAtom, jj's old neighbors can no longer be deleted
                    numNewDeleteCandidates-=numDeleteCandidateNbrs[jj];
                }
                else if (numNeighbors[jj] < 12 && deleteCand) {
                    // testAtom would be a deleteCandidate due to jj
//                        deleteCandidateTimes[testAtom.getLeafIndex()]++;
                    numNewDeleteCandidates++;
                }                        
            }
        }
        nbrs = potentialMaster.getNeighborManager(box).getDownList(testAtom)[0];
        for (int j=0; j<nbrs.getAtomCount(); j++) {
            IAtom jAtom = nbrs.getAtom(j);
            int jj = jAtom.getLeafIndex();
            dr.Ev1Mv2(pi, jAtom.getPosition());
            box.getBoundary().nearestImage(dr);
            double r2 = dr.squared();
            if (r2 < maxDistance*maxDistance) {
                nTestNbrs++;
                boolean deleteCand = false;
                if (r2 < maxInsertNbrDistance*maxInsertNbrDistance && r2 > minInsertNbrDistance*minInsertNbrDistance) {
                    for (int k=0; k<nbrVectors.length; k++) {
                        double s2 = dr.Mv1Squared(nbrVectors[k]);
                        if (s2 < maxInsertDistance*maxInsertDistance) {
                            deleteCand = true;
                            break;
                        }
                    }
                    if (deleteCand) {
                        nTestNbrsDeletion++;
                    }
                }
                if (numNeighbors[jj] == 12) {
                    // by inserting testAtom, jj's old neighbors can no longer be deleted
                    numNewDeleteCandidates-=numDeleteCandidateNbrs[jj];
                }
                else if (numNeighbors[jj] < 12 && deleteCand) {
                    // jj now has 1 more neighbor
                    // testAtom would be a deleteCandidate due to jj
//                        deleteCandidateTimes[testAtom.getLeafIndex()]++;
                    numNewDeleteCandidates++;
                }                        
            }
        }
        
        // it could be that testAtom itself has 12 or fewer neighbors
        // if so, those neighbors would now be delete candidates
        if (nTestNbrs<=12) {
            numNewDeleteCandidates += nTestNbrsDeletion;
        }

        return true;
    }

    protected void findCandidates() {
        NeighborListManager nbrManager = potentialMaster.getNeighborManager(box);
        IBoundary boundary = box.getBoundary();
        int numAtoms = box.getLeafList().getAtomCount();
        if (numNeighbors.length < numAtoms) {
            numNeighbors = new int[numAtoms];
            numNeighborCandidatesOnDelete = new int[numAtoms];
            deleteCandidateNbrsOnDelete = new List[numAtoms];
            for (int i=0; i<numAtoms; i++) {
                deleteCandidateNbrsOnDelete[i] = new ArrayList<Integer>();
            }
            numDeleteCandidateNbrs = new int[numAtoms];
            deleteCandidateTimes = new int[numAtoms];
        }
        insertCandidates.clear();
        deleteCandidates.clear();
        totalDeleteCandidateTimes = 0;
        for (int i=0; i<numAtoms; i++) {
            deleteCandidateNbrsOnDelete[i].clear();
            numNeighbors[i] = numNeighborCandidatesOnDelete[i] = numDeleteCandidateNbrs[i] = deleteCandidateTimes[i] = 0;
        }
        for (int i=0; i<numAtoms; i++) {
            IAtom iAtom = box.getLeafList().getAtom(i);
            IVector pi = iAtom.getPosition();
            IAtomList nbrsUp = nbrManager.getUpList(iAtom)[0];
            for (int j=0; j<nbrsUp.getAtomCount(); j++) {
                dr.Ev1Mv2(pi, nbrsUp.getAtom(j).getPosition());
                boundary.nearestImage(dr);
                double r2 = dr.squared();
                if (r2 < maxDistance*maxDistance) {
                    numNeighbors[i]++;
                    numNeighbors[nbrsUp.getAtom(j).getLeafIndex()]++;
                }
            }
        }
        double maxInsertNbrDistance = nbrDistance + maxInsertDistance;
        double minInsertNbrDistance = nbrDistance - maxInsertDistance;
        if (maxInsertNbrDistance > maxDistance) {
            throw new RuntimeException("oops");
        }
        for (int i=0; i<numAtoms; i++) {
            if (numNeighbors[i] < 14) {
                // the neighbors of i may be candidates for deletion.  after deleting
                // one of its neighbors, i would have <12 neighbors
                IAtom iAtom = box.getLeafList().getAtom(i);
                IVector pi = iAtom.getPosition();
                IAtomList nbrs = nbrManager.getUpList(iAtom)[0];
                for (int j=0; j<nbrs.getAtomCount(); j++) {
                    IAtom jAtom = nbrs.getAtom(j);
                    int jj = jAtom.getLeafIndex();
                    dr.Ev1Mv2(pi, jAtom.getPosition());
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    if (r2 < maxDistance*maxDistance) {
                        if (numNeighbors[i] == 12) {
                            // if we delete jj then i becomes an "insertCandidate"
                            numNeighborCandidatesOnDelete[jj]++;
                        }
                        else if (numNeighbors[i] == 13) {
                            // if we delete jj, i's neighbors (those within appropriate distance) will become delete candidates
                            deleteCandidateNbrsOnDelete[jj].add(i);
                        }
                        if (r2 < maxInsertNbrDistance*maxInsertNbrDistance && r2 > minInsertNbrDistance*minInsertNbrDistance) {
                            boolean success = false;
                            for (int k=0; k<nbrVectors.length; k++) {
                                double s2 = dr.Mv1Squared(nbrVectors[k]);
                                if (s2 < maxInsertDistance*maxInsertDistance) {
                                    success = true;
                                    break;
                                }
                            }
                            if (success) {
                                numDeleteCandidateNbrs[i]++;
                                if (numNeighbors[i] < 13) {
                                    // we need to know how many times jj shows up as a delete candidate
                                    deleteCandidateTimes[jj]++;
                                    totalDeleteCandidateTimes++;
                                    if (deleteCandidates.contains(jj)) continue;
                                    deleteCandidates.add(jj);
                                }
                            }
                        }
                    }
                }
                nbrs = nbrManager.getDownList(iAtom)[0];
                for (int j=0; j<nbrs.getAtomCount(); j++) {
                    IAtom jAtom = nbrs.getAtom(j);
                    int jj = jAtom.getLeafIndex();
                    dr.Ev1Mv2(pi, jAtom.getPosition());
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    if (r2 < maxDistance*maxDistance) {
                        if (numNeighbors[i] == 12) {
                            // if we delete jj then i becomes an "insertCandidate"
                            numNeighborCandidatesOnDelete[jj]++;
                        }
                        else if (numNeighbors[i] == 13) {
                            // if we delete jj, i's neighbors (those within appropriate distance) will become delete candidates
                            deleteCandidateNbrsOnDelete[jj].add(i);
                        }
                        if (r2 < maxInsertNbrDistance*maxInsertNbrDistance && r2 > minInsertNbrDistance*minInsertNbrDistance) {
                            boolean success = false;
                            for (int k=0; k<nbrVectors.length; k++) {
                                double s2 = dr.Mv1Squared(nbrVectors[k]);
                                if (s2 < maxInsertDistance*maxInsertDistance) {
                                    success = true;
                                    break;
                                }
                            }
                            if (success) {
                                numDeleteCandidateNbrs[i]++;
                                if (numNeighbors[i] < 13) {
                                    // we need to know how many times jj shows up as a delete candidate
                                    deleteCandidateTimes[jj]++;
                                    totalDeleteCandidateTimes++;
                                    if (deleteCandidates.contains(jj)) continue;
                                    deleteCandidates.add(jj);
                                }
                            }
                        }
                    }
                }
                if (numNeighbors[i] < 12) {
                    // we will attempt to insert next to i
                    insertCandidates.add(i);
                }
            }
        }
        dirty = false;
//        System.out.println(numAtoms+" "+insertCandidates.size()+" "+deleteCandidates.size());
        lastStepCount = integrator.getStepCount();
    }

    public double getA() {

        // our candidate for deletion was listed deleteCandidateTimes times
        return ((double)totalDeleteCandidateTimes)/(totalDeleteCandidateTimes + numNewDeleteCandidates);
    }

    // We have 2 atoms
    public AtomIterator affectedAtoms() {
        AtomArrayList list = (AtomArrayList)affectedAtomIterator.getList();
        list.clear();
        list.add(chosenAtom);
        list.add(testAtom);
        return affectedAtomIterator;
    }

    // Superclass methods have no idea what we're doing...
    public void acceptNotify() {
        box.removeMolecule(chosenAtom.getParentGroup());
        testAtom = chosenAtom;
        dirty = true;
    }

    // Superclass methods have no idea what we're doing...
    public void myRejectNotify() {
        box.removeMolecule(testAtom.getParentGroup());
    }
}

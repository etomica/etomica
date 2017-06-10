/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.integrator.Integrator;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.util.IEvent;
import etomica.util.IListener;

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
public class MCMoveInsertDeleteVacancy extends MCMoveInsertDeleteBiased implements IListener {

    protected final Vector dest;
    protected final Vector dr;
    protected Integrator integrator;
    protected long lastStepCount;
    protected boolean dirty;
    protected double minDistance, maxDistance, maxInsertDistance;
    protected List<Integer> insertCandidates, deleteCandidates;
    protected Set<Integer> deleteCandidateSet;
    protected int[] numNeighbors, numNeighborCandidatesOnDelete, deleteCandidateTimes, numDeleteCandidateNbrs;
    protected int totalDeleteCandidateTimes;
    protected PotentialMasterList potentialMaster;
    protected int numNewDeleteCandidates;
    protected int forced = 0;
    protected double oldLnA, oldB, newLnA;
    protected final Vector oldPosition;

    public MCMoveInsertDeleteVacancy(PotentialMaster potentialMaster,
                                     IRandom random, Space _space, Integrator integrator, double nbrDistance, int maxN, int maxVacancy) {
        super(potentialMaster, random, _space, maxN-maxVacancy, maxN);
        this.potentialMaster = (PotentialMasterList)potentialMaster;
        dest = _space.makeVector();
        dr = _space.makeVector();
        this.integrator = integrator;
        minDistance = 0.95;
        maxDistance = nbrDistance;
        maxInsertDistance = 1.1;
        if (maxInsertDistance > maxDistance) {
            throw new RuntimeException("nbrDistance must be greater than maxInsert distance");
        }
        insertCandidates = new ArrayList<Integer>();
        deleteCandidates = new ArrayList<Integer>();
        deleteCandidateSet = new HashSet<Integer>();
        numNeighbors = new int[0];
        numNeighborCandidatesOnDelete = new int[0];
        lastStepCount = -1;
        dirty = true;
        lnbias = new double[0];
        oldPosition = _space.makeVector();
    }

    public void setMinDistance(double minDistance) {
        this.minDistance = minDistance;
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

    public boolean doTrial() {
//        energyMeter.setTarget((IAtom)null);
//        energyMeter.setTarget((IMolecule)null);
//        uOldFull = energyMeter.getDataAsScalar();
        if (integrator.getStepCount() > lastStepCount) {
            forced = 0;
        }
        int numAtoms = box.getLeafList().getAtomCount();
        if (forced==1) {
            insert = !insert;
            System.out.println("forcing "+(insert ? "insertion" : "deletion"));
            forced=2;
        }
        else {
            insert = (random.nextInt(2) == 0);
        }
//        insert = false;
        numNewDeleteCandidates = 0;
        if (dirty || lastStepCount < integrator.getStepCount()) findCandidates();
        if(insert) {
            if(!reservoir.isEmpty()) testMolecule = reservoir.remove(reservoir.getMoleculeCount()-1);
            else testMolecule = species.makeMolecule();
            IAtom testAtom = testMolecule.getChildList().getAtom(0);

            int nInsertCandidates = insertCandidates.size();
            if (nInsertCandidates == 0) return false;
            IAtom partner = box.getLeafList().getAtom(insertCandidates.get(random.nextInt(nInsertCandidates)));
//            System.out.println("inserting next to "+partner);
            uOld = 0;
            double r2 = 0;
            while (r2 < minDistance*minDistance) {
                dest.setRandomInSphere(random);
                dest.TE(maxInsertDistance);
                r2 = dest.squared();
            }
            testAtom.getPosition().E(partner.getPosition());
            testAtom.getPosition().PE(dest);
            if (forced==2) {
                testAtom.getPosition().E(oldPosition);
            }
            box.addMolecule(testMolecule);
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar();

            // inserting testAtom might change some deleteCandidates
            // some existing deleteCandidates with nbrs=12 might have 13 now
            // we also need to see how many times testAtom shows up as a neighbor
            // of a deleteCandidate
            deleteCandidateTimes[testAtom.getLeafIndex()] = 0;
            Vector pi = testAtom.getPosition();
            IAtomList nbrs = potentialMaster.getNeighborManager(box).getUpList(testAtom)[0];
            int nTestNbrs = 0, nTestNbrsDeletion = 0;
            for (int j=0; j<nbrs.getAtomCount(); j++) {
                IAtom jAtom = nbrs.getAtom(j);
                int jj = jAtom.getLeafIndex();
                dr.Ev1Mv2(pi, jAtom.getPosition());
                box.getBoundary().nearestImage(dr);
                r2 = dr.squared();
                if (r2 < maxDistance*maxDistance) {
                    nTestNbrs++;
                    if (r2 > minDistance*minDistance && r2 < maxInsertDistance*maxInsertDistance) {
                        nTestNbrsDeletion++;
                    }
                    if (numNeighbors[jj] == 12) {
                        // jj now has 13 neighbors
                        // by inserting testAtom, jj's old neighbors can no longer be deleted
                        numNewDeleteCandidates-=numDeleteCandidateNbrs[jj];
                    }
                    else if (numNeighbors[jj] < 12 && r2 > minDistance*minDistance && r2 < maxInsertDistance*maxInsertDistance) {
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
                r2 = dr.squared();
                if (r2 < maxDistance*maxDistance) {
                    nTestNbrs++;
                    if (r2 > minDistance*minDistance && r2 < maxInsertDistance*maxInsertDistance) {
                        nTestNbrsDeletion++;
                    }
                    if (numNeighbors[jj] == 12) {
                        // by inserting testAtom, jj's old neighbors can no longer be deleted
                        numNewDeleteCandidates-=numDeleteCandidateNbrs[jj];
                    }
                    else if (numNeighbors[jj] < 12 && r2 > minDistance*minDistance && r2 < maxInsertDistance*maxInsertDistance) {
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
            for (int i : deleteCandidates) {
                icount += deleteCandidateTimes[i];
                if (icount > irand) {
                    ip = i;
                    break;
                }
            }
            if (forced==2) {
                ip = numAtoms-1;
            }
//            System.out.println("deleting "+ip);
            IAtom testAtom = box.getLeafList().getAtom(ip);
            testMolecule = testAtom.getParentGroup();
            //delete molecule only upon accepting trial
            energyMeter.setTarget(testMolecule);
            uOld = energyMeter.getDataAsScalar();
            
            // by deleting testAtom, we might turn other atoms into "insertCandidates"
            // numNeighborCandidates is how many of our neighbors have 12 nbrs.
            uNew = 0;
        } 
        return true;
    }
    
    protected void findCandidates() {
        NeighborListManager nbrManager = potentialMaster.getNeighborManager(box);
        Boundary boundary = box.getBoundary();
        int numAtoms = box.getLeafList().getAtomCount();
        numNeighbors = new int[numAtoms];
        numNeighborCandidatesOnDelete = new int[numAtoms];
        numDeleteCandidateNbrs = new int[numAtoms];
        deleteCandidateTimes = new int[numAtoms+1];
        insertCandidates.clear();
        deleteCandidateSet.clear();
        deleteCandidates.clear();
        totalDeleteCandidateTimes = 0;
        for (int i=0; i<numAtoms; i++) {
            IAtom iAtom = box.getLeafList().getAtom(i);
            Vector pi = iAtom.getPosition();
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
        for (int i=0; i<numAtoms; i++) {
            if (numNeighbors[i] < 13) {
                // the neighbors of i may be candidates for deletion.  after deleting
                // one of its neighbors, i would have <12 neighbors
                IAtom iAtom = box.getLeafList().getAtom(i);
                Vector pi = iAtom.getPosition();
                IAtomList nbrs = nbrManager.getUpList(iAtom)[0];
                for (int j=0; j<nbrs.getAtomCount(); j++) {
                    IAtom jAtom = nbrs.getAtom(j);
                    int jj = jAtom.getLeafIndex();
                    dr.Ev1Mv2(pi, jAtom.getPosition());
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    if (numNeighbors[i] == 12 && r2 < maxDistance*maxDistance) {
                        // if we delete jj then i becomes an "insertCandidate"
                        numNeighborCandidatesOnDelete[jj]++;
//                        if (jj==218) {
//                            System.out.println(i+" will become an insertCandidate if "+jj+" is deleted => "+numNeighborCandidatesOnDelete[jj]);
//                        }
                    }
                    if (r2 < maxInsertDistance*maxInsertDistance && r2 > minDistance*minDistance) {
                        // we need to know how many times jj shows up as a delete candidate
                        deleteCandidateTimes[jj]++;
                        numDeleteCandidateNbrs[i]++;
                        totalDeleteCandidateTimes++;
                        if (deleteCandidateSet.contains(jj)) continue;
                        deleteCandidates.add(jj);
                        deleteCandidateSet.add(jj);
                    }
                }
                nbrs = nbrManager.getDownList(iAtom)[0];
                for (int j=0; j<nbrs.getAtomCount(); j++) {
                    IAtom jAtom = nbrs.getAtom(j);
                    int jj = jAtom.getLeafIndex();
                    dr.Ev1Mv2(pi, jAtom.getPosition());
                    boundary.nearestImage(dr);
                    double r2 = dr.squared();
                    if (numNeighbors[i] == 12 && r2 < maxDistance*maxDistance) {
                        // if we delete jj then i becomes an "insertCandidate"
                        numNeighborCandidatesOnDelete[jj]++;
//                        if (jj==218) {
//                            System.out.println(i+" will become an insertCandidate if "+jj+" is deleted => "+numNeighborCandidatesOnDelete[jj]);
//                        }
                    }
                    if (r2 < maxInsertDistance*maxInsertDistance && r2 > minDistance*minDistance) {
                        // we need to know how many times jj shows up as a delete candidate
                        deleteCandidateTimes[jj]++;
                        numDeleteCandidateNbrs[i]++;
                        totalDeleteCandidateTimes++;
                        if (deleteCandidateSet.contains(jj)) continue;
                        deleteCandidates.add(jj);
                        deleteCandidateSet.add(jj);
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
        double lna = getLnBiasDiff();

        double shellV = 4.0/3.0*Math.PI*(maxInsertDistance*maxInsertDistance*maxInsertDistance - minDistance*minDistance*minDistance);
        double c = 0;
        if (insert) {
            //c = insertCandidates.size()*shellV*deleteCandidateTimes[testMolecule.getIndex()]/(totalDeleteCandidateTimes + numNewDeleteCandidates);
            c = insertCandidates.size()*shellV/(totalDeleteCandidateTimes + numNewDeleteCandidates);
//            System.out.println(insert+" c "+shellV+" "+totalDeleteCandidateTimes+" "+numNewDeleteCandidates+" "+insertCandidates.size());
//            System.out.print("insert candidates: ");
//            for (int i : insertCandidates) {
//                System.out.print(" "+i);
//            }
//            System.out.print("\n");
        }
        else {
            // our candidate for deletion was listed deleteCandidateTimes times
            int newInsertCandidates = insertCandidates.size() + numNeighborCandidatesOnDelete[testMolecule.getIndex()];
            if (insertCandidates.contains(testMolecule.getIndex())) newInsertCandidates--;
            c = totalDeleteCandidateTimes/(shellV*newInsertCandidates);
//            System.out.println(insert+" c "+shellV+" "+totalDeleteCandidateTimes+" "+newInsertCandidates+" "+insertCandidates.size()+" "+numNeighborCandidatesOnDelete[testMolecule.getIndex()]);
//            System.out.print("insert candidates: ");
//            for (int i : insertCandidates) {
//                System.out.print(" "+i);
//            }
//            System.out.print("\n");
        }
        if (forced == 2) {
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
//        System.out.println(insert+" log(c) = "+lc);
        return Math.exp(lna)*c;
    }

    public double getB() {
        double b = uOld - uNew;
        if (false) System.out.println(insert+" b = "+b);
        if (forced==2) {
            System.out.println("**** b = "+b);
            System.out.println("*old b = "+oldB);
            if (Math.abs((b+oldB)/(1+b-oldB)) > 1e-6) {
                throw new RuntimeException("oops");
            }
            forced = 0;
            return -Double.POSITIVE_INFINITY;
        }
        oldB = b;
        return b;
    }

    public void myAcceptNotify() {
//        System.out.println(insert+" "+Math.log(getA())+" "+getB()/0.7+" "+uOld+" "+uNew);
//        System.out.println("accepted "+(insert ? "insertion" : "deletion"));
        if (!insert) {
            oldPosition.E(testMolecule.getChildList().getAtom(0).getPosition());
        }
        super.myAcceptNotify();
        dirty = true;
//        System.out.println("accepted => "+box.getLeafList().getAtomCount());
//        forced = 1;
    }
    
    public void actionPerformed(IEvent event) {
        if (event instanceof MCMoveTrialCompletedEvent && ((MCMoveEvent)event).getMCMove() != this && ((MCMoveTrialCompletedEvent)event).isAccepted()) {
            dirty = true;
        }
    }

}

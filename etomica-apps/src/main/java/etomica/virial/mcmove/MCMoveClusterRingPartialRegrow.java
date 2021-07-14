/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * MCMove that partially regrows the beads of a ring polymer, accepting or
 * rejecting the move based on the sampling weight.  The move can (optionally)
 * regrow the beads such that the beads from multiple molecules are combined
 * to form a larger ring.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterRingPartialRegrow extends MCMoveBox {

    public MCMoveClusterRingPartialRegrow(PotentialMaster potentialMaster, IRandom random, Space _space) {
        this(potentialMaster, random, _space, new int[0][0]);
    }
    
    public MCMoveClusterRingPartialRegrow(PotentialMaster potentialMaster, IRandom random, Space _space, int[][] tangledMolecules) {
        super(potentialMaster);
        this.space = _space;
        this.random = random;
        this.tangledMolecules = tangledMolecules;
        setNumTrial(10);
        dcom = space.makeVector();
        leafIterator = new AtomIteratorLeafAtoms();
        myAtoms = new AtomArrayList();
        energyMeter = new MeterPotentialEnergy(potential);
	}
    
    public void setNumTrial(int newNumTrial) {
        nTrial = newNumTrial;
        rTrial = new Vector[nTrial];
        for (int i=0; i<nTrial; i++) {
            rTrial[i] = space.makeVector();
        }
        pkl = new double[nTrial];
    }

    public int getNumTrial() {
        return nTrial;
    }
    
    public void setEnergyFactor(double factor) {
        fac = factor;
    }
    
    public void setNumBeads(int newNumBeads) {
        maxNumBeads = newNumBeads;
        oldPositions = new Vector[maxNumBeads];
        for (int i=0; i<maxNumBeads; i++) {
            oldPositions[i] = space.makeVector();
        }
    }
    
    public int getNumBeads() {
        return numBeads;
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        int nMolecules = box.getMoleculeList().size();
        hist = new HistogramExpanding[nMolecules][0];
        for (int i=0; i<nMolecules; i++) {
            int nAtoms = box.getMoleculeList().get(i).getChildList().size();
            hist[i] = new HistogramExpanding[nAtoms];
            for (int j=0; j<nAtoms; j++) {
                hist[i][j] = new HistogramExpanding(0.04);
            }
        }
        leafIterator.setBox(p);
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
		uOld = energyMeter.getDataAsScalar();
		
        weightOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList molecules = box.getMoleculeList();
        wNew = 1;
        wOld = 1;
        double sigma = Math.sqrt(0.5/fac);
        numBeads = Math.round(random.nextInt((int)Math.round(maxNumBeads*0.9))+Math.round(maxNumBeads*0.1));

            iMolecule = random.nextInt(molecules.size());
            int i = iMolecule;
            atoms = null;
            int nAtoms = 0;
            boolean single = true;
            int[] tangled = null;
            for (int j=0; single && j<tangledMolecules.length; j++) {
                for (int k=0; k<tangledMolecules[j].length; k++) {
                    if (i==tangledMolecules[j][k]) {
                        single = false;
                        tangled = tangledMolecules[j];
                        i = tangled[0];
                        break;
                    }
                }
            }
            if (single) {
                atoms = molecules.get(i).getChildList();
                nAtoms = atoms.size();
            }
            else {
                myAtoms.clear();
                for (int j=0; j<tangled.length; j++) {
                    IAtomList jAtoms = molecules.get(tangled[j]).getChildList();
                    myAtoms.addAll(jAtoms);
                    nAtoms += jAtoms.size();
                }
                atoms = myAtoms;
            }
            dcom.E(0);

            kStart = random.nextInt(nAtoms);

            IAtom atom0 = atoms.get(kStart);
            Vector prevAtomPosition = atom0.getPosition();
            int kEnd = (kStart + numBeads + 1) % nAtoms;
            IAtom atomN = atoms.get(kEnd);
            Vector lastAtomPosition = atomN.getPosition();

            double pPrev = 1;
            int k = kStart;
            if (false) {
                // configurational bias consideration of the reverse move
                for (int m=0; m<numBeads; m++) {
                    k++;
                    if (k == nAtoms) {
                        k = 0;
                    }
                    IAtom kAtom = atoms.get(k);
                    double pSum = 0;
                    for (int l=0; l<nTrial-1; l++) {
                        for (int j=0; j<3; j++) {
                            rTrial[l].setX(j, sigma*random.nextGaussian());
                        }
                        rTrial[l].PE(prevAtomPosition);
                        double rin2 = rTrial[l].Mv1Squared(lastAtomPosition)/(numBeads-m);
                        pSum += Math.exp(-fac*rin2);
                    }
                    double rin2 = kAtom.getPosition().Mv1Squared(lastAtomPosition)/(numBeads-m);
                    double pOld = Math.exp(-fac*rin2);
                    wOld *= (pSum+pOld)/pPrev/nTrial;
                    pPrev = pOld;
                    prevAtomPosition = kAtom.getPosition();
                }
            }
            
            k = kStart;
            pPrev = 1;
            prevAtomPosition = atom0.getPosition();
            for (int m=0; m<numBeads; m++) {
                k++;
                if (k == nAtoms) {
                    k = 0;
                }
                IAtom kAtom = atoms.get(k);
                Vector kPosition = kAtom.getPosition();
                dcom.ME(kPosition);
                oldPositions[m].E(kPosition);

                double k1 = fac/(numBeads-m);
                double k2 = fac;
                double x0 = 1.0/(1.0 + 1.0/(numBeads-m));
                sigma = Math.sqrt(0.5/(k1 + k2));
                
                for (int j=0; j<3; j++) {
                    kPosition.setX(j, sigma*random.nextGaussian());
                }
                
                kPosition.PEa1Tv1(x0, prevAtomPosition);
                kPosition.PEa1Tv1(1.0-x0, lastAtomPosition);
                
                if (false ) {
                    // configurational bias consideration of the forward move
                    for (int l=0; l<nTrial; l++) {
                        for (int j=0; j<3; j++) {
                            rTrial[l].setX(j, sigma*random.nextGaussian());
                        }
                        rTrial[l].PE(prevAtomPosition);
                        double rin2 = rTrial[l].Mv1Squared(lastAtomPosition)/(numBeads-m);
                        pkl[l] = Math.exp(-fac*rin2); //*rin2;
    //                    if (k==1) System.out.println(rTrial[l].squared()+" "+(nAtoms-k)+" "+pkl[l]);
                        if (l>0) pkl[l] += pkl[l-1];
                    }
                    double ranl = random.nextDouble() * pkl[nTrial-1];
                    int chosenTrial = nTrial-1;
                    for (int l=0; l<nTrial-1; l++) {
                        if (pkl[l] >= ranl) {
                            chosenTrial = l;
                            break;
                        }
                    }
    //                if (k==1) System.out.println(rTrial[chosenTrial].squared());
                    kAtom.getPosition().E(rTrial[chosenTrial]);
                    wNew *= pkl[nTrial-1]/pPrev/nTrial;
                    pPrev = pkl[chosenTrial] - (chosenTrial == 0 ? 0 : pkl[chosenTrial-1]);
                }
                dcom.PE(kPosition);
                prevAtomPosition = kPosition;
            }
            dcom.TE(-1.0/nAtoms);

            for (k=0; k<nAtoms; k++) {
                atoms.get(k).getPosition().PE(dcom);
            }

        ((BoxCluster)box).trialNotify();
        weightNew = calcWeight();
//        System.out.println(wOld+" =?=> "+wNew);
        
        
        uNew = energyMeter.getDataAsScalar();
        
		return true;
	}
	
	protected double calcWeight() {
		return ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
	}

    public double getChi(double temperature) {
        return wNew / wOld * weightNew / weightOld * Math.exp(-(uNew - uOld) / temperature);
    }
    
    public void rejectNotify() {
        int nAtoms = atoms.size();
        for (int k=0; k<nAtoms; k++) {
            atoms.get(k).getPosition().ME(dcom);
        }
        for (int k=kStart, m=0; m<numBeads; m++) {
            k++;
            if (k == nAtoms) {
                k = 0;
            }
            IAtom kAtom = atoms.get(k);
            kAtom.getPosition().E(oldPositions[m]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
//        System.out.println("accepted");
    	((BoxCluster)box).acceptNotify();
    }
    
    public double energyChange() {
        return 0;
    }
    
    public AtomIterator affectedAtoms() {
        return leafIterator;
    }
    
    private static final long serialVersionUID = 1L;
    protected IAtomList atoms;
    protected final Space space;
    protected final IRandom random;
    protected Vector[] oldPositions;
    protected Vector[] rTrial;
    protected int nTrial;
    protected double[] pkl;
    // Rosenbluth weights
    protected double wOld, wNew;
    // cluster weights
    protected double weightOld, weightNew, uOld, uNew;
    protected final Vector dcom;
    protected final AtomIteratorLeafAtoms leafIterator;
    protected double fac;
    public HistogramExpanding[][] hist;
    protected final int[][] tangledMolecules;
    protected final AtomArrayList myAtoms;
    protected int numBeads, maxNumBeads;
    protected int kStart, iMolecule;
    protected final MeterPotentialEnergy energyMeter;
}

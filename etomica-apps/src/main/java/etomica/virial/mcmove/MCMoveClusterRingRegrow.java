/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * MCMove that fully regrows the beads of a ring polymer, accepting or
 * rejecting the move based on the sampling weight.  The move can (optionally)
 * regrow the beads such that the beads from multiple molecules are combined
 * to form a larger ring.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterRingRegrow extends MCMoveBox {

    public MCMoveClusterRingRegrow(IRandom random, Space _space) {
        this(random, _space, new int[0][0]);
    }
    
    public MCMoveClusterRingRegrow(IRandom random, Space _space, int[][] tangledMolecules) {
        super();
        this.space = _space;
        this.random = random;
        this.tangledMolecules = tangledMolecules;
        setNumTrial(10);
        com = space.makeVector();
        com0 = space.makeVector();
        myAtoms = new AtomArrayList();
        wOld = 1e-10;
	}
    
    /**
     * Sets the number of configurational bias trials
     */
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
    
    /**
     * sets the harmonic bond "energy" factor
     * The probability of the bond have length x is proportional to exp(-factor*x^2)
     * (there is no temperature involved here)
     */
    public void setEnergyFactor(double factor) {
        fac = factor;
    }

    public void setBox(Box p) {
        super.setBox(p);
        int nMolecules = box.getMoleculeList().size();
        oldPositions = new Vector[nMolecules][0];
        for (int i=0; i<nMolecules; i++) {
            int nAtoms = box.getMoleculeList().get(i).getChildList().size();
            oldPositions[i] = new Vector[nAtoms];
            for (int j=0; j<nAtoms; j++) {
                oldPositions[i][j] = space.makeVector();
            }
        }
    }

	public boolean doTrial() {
        weightOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList molecules = box.getMoleculeList();
        // Rosenbluth weight
        wNew = 1;

        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = null;
            int nAtoms = 0;
            boolean single = true;
            boolean skip = false;
            int[] tangled = null;
            // this is to handle the "exchange" case where chains from two or
            // more molecules are connected together to form a larger ring.
            // usually, we have to resort to partial regrowth
            for (int j=0; single && j<tangledMolecules.length; j++) {
                for (int k=0; k<tangledMolecules[j].length; k++) {
                    if (i==tangledMolecules[j][k]) {
                        single = false;
                        skip = k>0;
                        tangled = tangledMolecules[j];
                        break;
                    }
                }
            }
            // skip this molecule (for now) if it is not the first in the list
            // of "tangled" molecules
            if (skip) continue;
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
            // determine the old center of mass
            com.E(0);
            for (int k=0; k<nAtoms; k++) {
                com.PE(atoms.get(k).getPosition());
            }
            com.TE(1.0/nAtoms);

            IAtom atom0 = atoms.get(0);
            oldPositions[i][0].E(atom0.getPosition());
            // put the first atom at the origin, as a starting point
            // we'll translate everything back to the original center of mass later
            atom0.getPosition().E(0);
            Vector prevAtomPosition = atom0.getPosition();
            double pPrev = 1;
            com0.E(0);
            IMolecule moleculei = molecules.get(i);
            int iTangled = 0;
            int kStart = 0;
            for (int k=1; k<nAtoms; k++) {
                if (k-kStart >= moleculei.getChildList().size()) {
                    kStart += moleculei.getChildList().size();
                    iTangled++;
                    i = tangled[iTangled];
                    moleculei = molecules.get(i);
                }
                IAtom kAtom = atoms.get(k);
                Vector kPosition = kAtom.getPosition();
                oldPositions[i][k-kStart].E(kPosition);

                double k1 = fac/(nAtoms-k);
                double k2 = fac;
                double K = k1 + k2;
                double sigma = Math.sqrt(0.5/K);
                for (int j=0; j<3; j++) {
                    kPosition.setX(j, sigma*random.nextGaussian());
                }
                
                kPosition.PEa1Tv1(k2/K, prevAtomPosition);

                if (false) {
                    // this handles configurational bias, which is (currently) unnecessary.
                    for (int l=0; l<nTrial; l++) {
                        for (int j=0; j<3; j++) {
                            rTrial[l].setX(j, sigma*random.nextGaussian());
                        }
                        rTrial[l].PE(prevAtomPosition);
                        double rin2 = rTrial[l].squared()/(nAtoms-k);
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
                    com0.PE(kAtom.getPosition());
                    wNew *= pkl[nTrial-1]/pPrev/nTrial;
                    pPrev = pkl[chosenTrial] - (chosenTrial == 0 ? 0 : pkl[chosenTrial-1]);
                    if (Double.isNaN(wNew)) {
                        throw new RuntimeException("oops");
                    }
                }
                com0.PE(kPosition);
                prevAtomPosition = kPosition;
            }
            com0.TE(1.0/nAtoms);
            com.ME(com0);

            // translate back to original COM
            for (int k=0; k<nAtoms; k++) {
                atoms.get(k).getPosition().PE(com);
            }
        }
        
		((BoxCluster)box).trialNotify();
        weightNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}

    public double getChi(double temperature) {
        // we skip calculation of wOld because we're the only intramolecular move in town.
        return wNew/wOld * weightNew/weightOld;
    }

    public void rejectNotify() {
        IMoleculeList molecules = box.getMoleculeList();

        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            int nAtoms = atoms.size();
            for (int j=0; j<nAtoms; j++) {
                atoms.get(j).getPosition().E(oldPositions[i][j]);
            }
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
        wOld = wNew;
    	((BoxCluster)box).acceptNotify();
    }
    
    public double energyChange() {
        return 0;
    }

    protected final Space space;
    protected final IRandom random;
    protected Vector[][] oldPositions;
    protected Vector[] rTrial;
    protected int nTrial;
    protected double[] pkl;
    // Rosenbluth weights
    protected double wOld, wNew;
    // cluster weights
    protected double weightOld, weightNew;
    protected final Vector com, com0;
    protected double fac;
    protected final int[][] tangledMolecules;
    protected final AtomArrayList myAtoms;
}

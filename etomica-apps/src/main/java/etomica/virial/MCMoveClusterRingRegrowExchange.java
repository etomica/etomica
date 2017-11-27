/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomHydrogen;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientation3D;
import etomica.util.random.IRandom;

/**
 * MCMove that fully regrows the beads of a ring polymer, accepting or
 * rejecting the move based on the sampling weight.  The move can (optionally)
 * regrow the beads such that the beads from multiple molecules are combined
 * to form a larger ring.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterRingRegrowExchange extends MCMoveBox {    
    
    public MCMoveClusterRingRegrowExchange(IRandom random, Space _space) {
        super(null);
        this.space = _space;
        this.random = random;        
        setNumTrial(10);
        com = space.makeVector();
        com0 = space.makeVector();
        leafIterator = new AtomIteratorLeafAtoms();
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
     * The probability of the bond to have length x is proportional to exp(-factor*x^2)
     * (there is no temperature involved here)
     */
    public void setEnergyFactor(double factor) {
        fac = factor;
    }

    public void setBox(Box p) {
        super.setBox(p);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        oldPositions = new Vector[nMolecules][0];
        oldOrientations = new IOrientation3D [nMolecules][0];
        oldBondlengths = new double [nMolecules][0];
        for (int i=0; i<nMolecules; i++) {
            int nAtoms = box.getMoleculeList().getMolecule(i).getChildList().getAtomCount();
            P = 2*nAtoms;
            oldOrientations[i] = new IOrientation3D[nAtoms];
            oldBondlengths[i] = new double [nAtoms];
            oldPositions[i] = new Vector[nAtoms];
            for (int j=0; j<nAtoms; j++) {
            	oldOrientations[i][j] = (IOrientation3D)space.makeOrientation();
                oldPositions[i][j] = space.makeVector();
            }
        }
        leafIterator.setBox(p);
    }

	public boolean doTrial() {
		if (P <= 0) throw new RuntimeException("value of P is not set");
        weightOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList molecules = box.getMoleculeList();
        // Rosenbluth weight
        wNew = 1;

        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();            
            // determine the old center of mass
            int nAtoms = atoms.getAtomCount();
            com.E(0);
            for (int k=0; k<nAtoms; k++) {
            	AtomHydrogen kAtom = ((AtomHydrogen)atoms.getAtom(k));                
                oldPositions[i][k].E(kAtom.getPosition());
                com.PE(oldPositions[i][k]);
                oldOrientations[i][k].setDirection(kAtom.getOrientation().getDirection());
                oldBondlengths[i][k] = kAtom.getBondLength();
            }
            com.TE(1.0/nAtoms);
            
            AtomHydrogen atom0 = (AtomHydrogen) atoms.getAtom(0);
            // put the first atom at the origin, as a starting point
            // we'll translate everything back to the original center of mass later
            atom0.getPosition().E(0);
            Vector prevAtomPosition = atom0.getPosition();
                        
            Vector[] newPositions = new Vector[P];
            newPositions[0] = space.makeVector();
            newPositions[0].E(0);
            for (int k=1; k<P; k++) {
            	newPositions[k] = space.makeVector();                
                double k1 = fac/(P-k);
                double k2 = fac;
                double K = k1 + k2;
                double sigma = Math.sqrt(0.5/K);
                for (int j=0; j<3; j++) {
                    newPositions[k].setX(j, sigma*random.nextGaussian());
                }                
                newPositions[k].PEa1Tv1(k2/K, prevAtomPosition);                
                com0.PE(newPositions[k]);
                prevAtomPosition = newPositions[k];
            }
            com0.TE(1.0/P);
            com.ME(com0);
            
            for (int j=0; j<nAtoms; j++) {//setting positions, bond lengths and orientations
            	AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);            	
            	Vector dummy = space.makeVector();
            	dummy.Ev1Pv2(newPositions[j], newPositions[j+nAtoms]);
            	dummy.TE(0.5);            	
            	jAtom.getPosition().E(dummy);
            	double bl = newPositions[j].Mv1Squared(newPositions[j+nAtoms]);
            	bl = Math.sqrt(bl);
            	jAtom.setBondLength(bl);
            	dummy.Ev1Mv2(newPositions[j], newPositions[j+nAtoms]);
            	dummy.normalize();
            	if(dummy.isNaN() || dummy.isZero()) throw new RuntimeException("oops!");
            	jAtom.getOrientation().setDirection(dummy);
            }

            // translate back to original COM
            for (int k=0; k<nAtoms; k++) {
                atoms.getAtom(k).getPosition().PE(com);
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

        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList atoms = molecules.getMolecule(i).getChildList();
            int nAtoms = atoms.getAtomCount();
            for (int j=0; j<nAtoms; j++) {
                atoms.getAtom(j).getPosition().E(oldPositions[i][j]);
                ((AtomHydrogen)atoms.getAtom(j)).setBondLength(oldBondlengths[i][j]);
                ((AtomHydrogen)atoms.getAtom(j)).getOrientation().E(oldOrientations[i][j]);;
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
    
    public AtomIterator affectedAtoms() {
        return leafIterator;
    }
    
    private static final long serialVersionUID = 1L;
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
    protected final AtomIteratorLeafAtoms leafIterator;
    protected double fac;    
    protected final AtomArrayList myAtoms;
    protected int P = -1;
    protected IOrientation3D[][] oldOrientations;
    protected double[][] oldBondlengths;
}

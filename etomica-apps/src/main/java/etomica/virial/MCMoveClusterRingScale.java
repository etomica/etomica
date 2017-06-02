/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IAtomList;
import etomica.box.Box;
import etomica.api.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.space.Vector;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.Space;

/**
 * MCMove that scales the atoms of a molecule in or out.  This is helpful
 * for PIMC (XC) when the beads get stuck in a non-overlapping configuration.
 */
public class MCMoveClusterRingScale extends MCMoveBox {

    public MCMoveClusterRingScale(PotentialMaster potentialMaster, IRandom random, Space _space) {
        this(potentialMaster, random, _space, new int[0][0]);
    }
    
    public MCMoveClusterRingScale(PotentialMaster potentialMaster, IRandom random, Space _space, int[][] tangledMolecules) {
        super(potentialMaster);
        this.space = _space;
        this.random = random;
        this.tangledMolecules = tangledMolecules;
        com = space.makeVector();
        leafIterator = new AtomIteratorLeafAtoms();
        myAtoms = new AtomArrayList();
        energyMeter = new MeterPotentialEnergy(potential);
	}
    
    public void setEnergyFactor(double factor) {
        fac = factor;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        leafIterator.setBox(p);
        energyMeter.setBox(p);
    }
    
	public boolean doTrial() {
		uOld = energyMeter.getDataAsScalar();
		
        weightOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList molecules = box.getMoleculeList();

        iMolecule = random.nextInt(molecules.getMoleculeCount());
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
            atoms = molecules.getMolecule(i).getChildList();
            nAtoms = atoms.getAtomCount();
        }
        else {
            myAtoms.clear();
            for (int j=0; j<tangled.length; j++) {
                IAtomList jAtoms = molecules.getMolecule(tangled[j]).getChildList();
                myAtoms.addAll(jAtoms);
                nAtoms += jAtoms.getAtomCount();
            }
            atoms = myAtoms;
        }
        com.E(0);
        for (int j=0; j<nAtoms; j++) {
            com.PE(atoms.getAtom(j).getPosition());
        }
        com.TE(1.0/nAtoms);
        //System.out.println(com);
        boolean scaleUp = random.nextInt(2) == 1;
        scale = scaleUp ? 1.001 : 1.0/1.001;
        double uOld = 0;
        for (int j=0; j<nAtoms; j++) {
            Vector p = atoms.getAtom(j).getPosition();
            double uj = 0;
            if (j==0) {
                uj += p.Mv1Squared(atoms.getAtom(nAtoms-1).getPosition());
            }
            if (j<nAtoms-1) {
                uj += p.Mv1Squared(atoms.getAtom(j+1).getPosition());
            }
//            System.out.println(j+" "+uj*fac+" "+uOld*fac);
            uOld += uj;
            p.ME(com);
            p.TE(scale);
            p.PE(com);
        }
        double uNew = scale*scale*uOld;
        wRatio = Math.exp(-fac*(uNew-uOld) + 3*(nAtoms-1)*Math.log(scale));

        ((BoxCluster)box).trialNotify();
        weightNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
//        System.out.println(wOld+" =?=> "+wNew);
        
        uNew = energyMeter.getDataAsScalar();
        
		return true;
	}
	
    public double getA() {
        // we skip calculation of wOld because we're the only intramolecular move in town.
        return wRatio * weightNew/weightOld;
    }

    public double getB() {
        return -(uNew - uOld);
    }
    
    public void rejectNotify() {
        int nAtoms = atoms.getAtomCount();
        for (int j=0; j<nAtoms; j++) {
            Vector p = atoms.getAtom(j).getPosition();
            double uj = 0;
            if (j==0) {
                uj += p.Mv1Squared(atoms.getAtom(nAtoms-1).getPosition());
            }
            if (j<nAtoms-1) {
                uj += p.Mv1Squared(atoms.getAtom(j+1).getPosition());
            }
            p.ME(com);
            p.TE(1.0/scale);
            p.PE(com);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }
    
    public double energyChange() {
        return 0;
    }
    
    public AtomIterator affectedAtoms() {
        return leafIterator;
    }
    
    private static final long serialVersionUID = 1L;
    protected double scale;
    protected double wRatio;
    protected IAtomList atoms;
    protected final Space space;
    protected final IRandom random;
    protected double weightOld, weightNew, uOld, uNew;
    protected final Vector com;
    protected final AtomIteratorLeafAtoms leafIterator;
    protected double fac;
    protected final int[][] tangledMolecules;
    protected final AtomArrayList myAtoms;
    protected int iMolecule;
    protected final MeterPotentialEnergy energyMeter;
}

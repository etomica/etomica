/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.paralleltempering;

import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPT;
import etomica.integrator.mcmove.MCMove;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.virial.BoxCluster;

/**
 * Swaps configurations and pairSet between boxes for a virial clustering simulation. 
 */
public class MCMoveSwapCluster extends MCMove implements IntegratorPT.MCMoveSwap {

    private static final long serialVersionUID = 1L;
    private IntegratorMC integrator1, integrator2;
    private AtomIteratorLeafAtoms iterator1 = new AtomIteratorLeafAtoms();
    private AtomIteratorLeafAtoms iterator2 = new AtomIteratorLeafAtoms();
    private Vector r;
    private BoxCluster box1, box2;
    private double weightOld1, weightOld2;
    private double weightNew1, weightNew2;
    private final Box[] swappedBoxes = new Box[2];

    public MCMoveSwapCluster(IntegratorMC integrator1, IntegratorMC integrator2, Space space) {
        super(null);
        r = space.makeVector();
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;		
    }

    public boolean doTrial() {
        if(box1 == null || box2 == null) {
            box1 = (BoxCluster)integrator1.getBox();
            box2 = (BoxCluster)integrator2.getBox();
            iterator1.setBox(box1);
            iterator2.setBox(box2);
        }

        weightOld1 = box1.getSampleCluster().value(box1);
        weightOld2 = box2.getSampleCluster().value(box2);
        
//        System.out.println("in trial "+integrator2.getWeight()+" "+weightOld2);
//        System.out.println("in trial "+integrator1.getWeight()+" "+weightOld1);
        iterator1.reset();
        iterator2.reset();

        for (IAtom a1 = iterator1.nextAtom(); a1 != null;
             a1 = iterator1.nextAtom()) {
            IAtom a2 = iterator2.nextAtom();

            //swap coordinates
            r.E(a1.getPosition());
            
            a1.getPosition().E(a2.getPosition());
            a2.getPosition().E(r);
        }

        //assumes energy will be determined using only pairSets in boxes
        box1.trialNotify();
        box2.trialNotify();
		
        weightNew1 = weightNew2 = Double.NaN;
        return true;
    }

    public double getChi(double temperature) {
        weightNew1 = box1.getSampleCluster().value(box1);
        weightNew2 = box2.getSampleCluster().value(box2);
//        System.out.println(weightOld1+" "+weightOld2+" "+weightNew1+" "+weightNew2);
        return  (weightNew1 * weightNew2) / (weightOld1 * weightOld2);
    }
	
    /**
     * Swaps positions of molecules in two boxes.
     */
    public void acceptNotify() {
//        System.out.println("accepted");
		
        box1.acceptNotify();
        box2.acceptNotify();
    }
	
    public void rejectNotify() {
//        System.out.println("rejected");
        iterator1.reset();
        iterator2.reset();

        for (IAtom a1 = iterator1.nextAtom(); a1 != null;
             a1 = iterator1.nextAtom()) {
            IAtom a2 = iterator2.nextAtom();

            //swap coordinates
            r.E(a1.getPosition());
            
            a1.getPosition().E(a2.getPosition());
            a2.getPosition().E(r);
        }

        box1.rejectNotify();
        box2.rejectNotify();
    }
    
    public double energyChange(Box box) {
        if(box == box1) return weightNew1/weightOld1;
        if(box == box2) return weightNew2/weightOld2;
        return 0.0;
    }

    /**
     * Implementation of MCMoveSwap interface
     */
    public Box[] swappedBoxes() {
        swappedBoxes[0] = box1;
        swappedBoxes[1] = box2;
        return swappedBoxes;
    }

    public AtomIterator affectedAtoms(Box p) {
        if(p == box1) {
            return iterator1;
        }
        else if (p == box2) {
            return iterator2;
        }
        return AtomIteratorNull.INSTANCE;
    }
	
    public final static SwapFactory FACTORY = new SwapFactory();
    
    protected static class SwapFactory implements IntegratorPT.MCMoveSwapFactory, java.io.Serializable {
        public MCMove makeMCMoveSwap(IntegratorBox integrator1, IntegratorBox integrator2, Space _space) {
            return new MCMoveSwapCluster((IntegratorMC)integrator1, (IntegratorMC)integrator2, _space);
        }
        private static final long serialVersionUID = 1L;
    } 
	
}


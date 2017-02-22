/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.entropylottery;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.ISpace;


/**
 * Monte Carlo move that moves an atom by +/- 1 unit in a random dimension.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomAdjacent extends MCMoveBox {
    
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected IVectorMutable translationVector;
    protected IAtom atom;
    protected AtomSource atomSource;
    protected final IRandom random;
    private final ISpace space;

    public MCMoveAtomAdjacent(IRandom random, ISpace _space) {
        super(null);
        this.random = random;
        this.space = _space;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        perParticleFrequency = true;
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        translationVector.E(0);
        int i = random.nextInt(translationVector.getD());
        translationVector.setX(i, random.nextInt(2)*2-1);
        atom.getPosition().PE(translationVector);
        return true;
    }
    
    public double getA() {return 1.0;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        IVectorMutable position = atom.getPosition();
        IVector dimensions = box.getBoundary().getBoxSize();
        for (int i=0; i<position.getD(); i++) {
            // if we're non-periodic, ensure we didn't try to jump over the boundary
            int x = (int)Math.round(position.getX(i)+dimensions.getX(i)*0.5-0.5);
            if (x < 0 || x >= (int)Math.round(dimensions.getX(i))) {
                if (!box.getBoundary().getPeriodicity(i)) {
                    // failure
                    return Double.NEGATIVE_INFINITY;
                }
                //wrap around -- OK, it's a hack.  deal with it.
                if (x < 0) {
                    position.setX(i, position.getX(i)+dimensions.getX(i));
                }
                else {
                    position.setX(i, 0);
                }
            }
        }
        return 0;
    }
    
    public double energyChange() {return 0;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
        translationVector.TE(-1);
        atom.getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        translationVector = space.makeVector();
        atomSource.setBox(p);
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }

    private static final long serialVersionUID = 1L;
}
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.ising;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 */
public class MCMoveSpinFlip extends MCMoveBox {

    private static final long serialVersionUID = 1L;
    protected final IRandom random;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;

    /**
     * @param potentialMaster
     * @param nBoxs
     */
    public MCMoveSpinFlip(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#doTrial()
     */
    public boolean doTrial() {
        IAtomList leafList = box.getLeafList();
        atom = leafList.get(random.nextInt(leafList.size()));
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        atom.getPosition().TE(-1);
        uNew = Double.NaN;
        return true;
    }

    public double getChi(double temperature) {
        uNew = energyMeter.getDataAsScalar();
        return Math.exp(-(uNew - uOld) / temperature);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#acceptNotify()
     */
    public void acceptNotify() {
        //nothing to do
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#rejectNotify()
     */
    public void rejectNotify() {
        atom.getPosition().TE(-1);
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#affectedAtoms(etomica.Box)
     */
    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    /* (non-Javadoc)
     * @see etomica.integrator.MCMove#energyChange(etomica.Box)
     */
    public double energyChange() {
        return uNew - uOld;
    }
}

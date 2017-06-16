/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.box.Box;
import etomica.math.function.IFunction;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.simulation.Simulation;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Pressure;
import etomica.math.function.Function;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolume extends MCMoveBoxStep {

    protected double pressure;
    protected MeterPotentialEnergy energyMeter;
    protected BoxInflate inflate;
    protected final int D;
    protected IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected IFunction vBias;

    protected double biasOld, uOld, hOld, vNew, vScale, hNew;
    protected double uNew = Double.NaN;
    protected double temperature;

    public MCMoveVolume(Simulation sim, PotentialMaster potentialMaster,
                        Space _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0);
    }

    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolume(PotentialMaster potentialMaster, IRandom random,
                        Space space, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.D = space.D();
        inflate = new BoxInflate(space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        vBias = new Function.Constant(1);
    }

    public void setInflater(BoxInflate newInflate) {
        inflate = newInflate;
        inflate.setBox(box);
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }

    public void setVolumeBias(IFunction vBias) {
        this.vBias = vBias;
    }

    /**
     * Sets the temperature so that chi can be compted in getA, avoiding
     * potential roundoff.  B and the exponential argument in A can both be
     * very large for a large system, but they tend to cancel.  Returned
     * separately, they might be 0 and infinity, causing big trouble.  If 
     * the move knows the temperature, all can be computed in getA without
     * risk of roundoff.
     */
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    public boolean doTrial() {
        double vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        biasOld = vBias.f(vOld);
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        inflate.setScale(rScale);
        //cells+neighbords get updated here
        inflate.actionPerformed();
        uNew = energyMeter.getDataAsScalar();
        hNew = uNew + pressure*vNew;
        return true;
    }//end of doTrial

    public double getA() {
        // N, not N+1 here because of the shell volume
        // D. S. Corti, Mol. Phys. 100, 1887 (2002).
        double biasNew = vBias.f(box.getBoundary().volume());
        if (temperature != 0) {
            return biasNew/biasOld*Math.exp(box.getMoleculeList().getMoleculeCount()*vScale - (hNew-hOld)/temperature);
        }
        return biasNew/biasOld*Math.exp(box.getMoleculeList().getMoleculeCount()*vScale);
    }

    public double getB() {
        if (temperature > 0) return 0;
        return -(hNew - hOld);
    }

    public void acceptNotify() {  /* do nothing */}

    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}

    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.action.BoxInflate;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.math.function.Function;
import etomica.math.function.IFunction;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Pressure;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolume extends MCMoveBoxStep {

    protected final IntegratorBox integrator;
    protected double pressure;
    protected BoxInflate inflate;
    protected final int D;
    protected final IRandom random;
    protected IFunction vBias;

    protected double biasOld, uOld, hOld, vNew, vScale, hNew;
    protected double uNew = Double.NaN;

    public MCMoveVolume(IntegratorBox integrator, IRandom random,
                        double pressure) {
        super();
        this.integrator = integrator;
        this.random = random;
        Space space = integrator.getBox().getSpace();
        this.D = space.D();
        inflate = new BoxInflate(space);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        vBias = new Function.Constant(1);
    }

    public void setInflater(BoxInflate newInflate) {
        inflate = newInflate;
        inflate.setBox(box);
    }

    public void setBox(Box p) {
        super.setBox(p);
        inflate.setBox(p);
    }

    public void setVolumeBias(IFunction vBias) {
        this.vBias = vBias;
    }

    public boolean doTrial() {
        double vOld = box.getBoundary().volume();
        uOld = integrator.getPotentialEnergy();
        hOld = uOld + pressure * vOld;
        biasOld = vBias.f(vOld);
        vScale = (2. * random.nextDouble() - 1.) * stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale / D);
        inflate.setScale(rScale);
        //cells+neighbords get updated here
        inflate.actionPerformed();
        PotentialCompute potentialCompute = integrator.getPotentialCompute();
        potentialCompute.init();
        uNew = potentialCompute.computeAll(false);
        hNew = uNew + pressure * vNew;
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        // N, not N+1 here because of the shell volume
        // D. S. Corti, Mol. Phys. 100, 1887 (2002).
        double biasNew = vBias.f(box.getBoundary().volume());
        return biasNew / biasOld * Math.exp(box.getMoleculeList().size() * vScale - (hNew - hOld) / temperature);
    }

    public void acceptNotify() {  /* do nothing */}

    public void rejectNotify() {
        inflate.undo();
        PotentialCompute potentialCompute = integrator.getPotentialCompute();
        potentialCompute.init();
        potentialCompute.computeAll(false);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void setPressure(double p) {
        pressure = p;
    }

    public final double getPressure() {
        return pressure;
    }

    public Dimension getPressureDimension() {
        return Pressure.DIMENSION;
    }
}

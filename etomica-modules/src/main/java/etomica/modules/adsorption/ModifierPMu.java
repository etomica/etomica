/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.integrator.IntegratorBox;
import etomica.modifier.Modifier;
import etomica.potential.P2HardGeneric;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

public class ModifierPMu implements Modifier {

    protected final EOSSW eos;
    protected final MyMCMove move;
    protected final P2HardGeneric p2;
    protected final MeterExcessAdsorbed meterCountAtoms;
    protected double oldValue = 0;

    protected IntegratorBox integrator;

    public ModifierPMu(P2HardGeneric p2, IntegratorBox integrator, EOSSW eos, MyMCMove move, MeterExcessAdsorbed meterCountAtoms) {
        this.p2 = p2;
        this.integrator = integrator;
        this.move = move;
        this.eos = eos;
        this.meterCountAtoms = meterCountAtoms;
        eos.setTemperature(integrator.getTemperature());
        reset();
    }

    public void setValue(double newValue) {
//        System.out.println("newValue "+newValue);
        double temperature = integrator.getTemperature();
        eos.setTemperature(temperature);
//        double sigma = p2.getCoreDiameter();
        double p = Math.pow(10, newValue); // * eos.pSat() / (sigma*sigma*sigma);
        double mu = eos.muForPressure(p);
        meterCountAtoms.setPressure(p);
//        System.out.println(temperature+" "+p+" "+mu);
        move.setMu(mu);
        oldValue = newValue;
    }

    public void reset() {
        eos.setSigma(p2.getCollisionDiameter(0));
        eos.setEpsilon(-p2.getEnergyForState(1));
        eos.setLambda(p2.getCollisionDiameter(1) / p2.getCollisionDiameter(0));
    }

    public double getValue() {
        // cheat.
        return oldValue;
    }

    public Dimension getDimension() {
        return Null.DIMENSION;
    }

    public String getLabel() {
        return "pressure";
    }
}

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.integrator.IntegratorBox;
import etomica.modifier.Modifier;
import etomica.units.dimensions.Dimension;
import etomica.units.Kelvin;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.virial.ClusterAbstract;
import etomica.virial.MCMoveClusterRingRegrow;

public class ModifierTPI implements Modifier {

    protected final MCMoveClusterRingRegrow ring0, ring1;
    protected final ClusterAbstract targetCluster;
    protected final double heMass;
    protected final int nBeads;
    protected double oldValue = 0;
    
    protected IntegratorBox integrator;
    
    public ModifierTPI(double heMass, int nBeads, MCMoveClusterRingRegrow ring0, MCMoveClusterRingRegrow ring1, ClusterAbstract targetCluster) {
        this.heMass = heMass;
        this.nBeads = nBeads;
        this.ring0 = ring0;
        this.ring1 = ring1;
        this.targetCluster = targetCluster;
    }
    
    public void setValue(double newValue) {
        double temperature = Kelvin.UNIT.toSim(Math.pow(10, newValue));
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*heMass*temperature);
        ring0.setEnergyFactor(nBeads*Math.PI/(lambda*lambda));
        ring1.setEnergyFactor(nBeads*Math.PI/(lambda*lambda));
        targetCluster.setTemperature(temperature);
        oldValue = newValue;
    }
    
    public double getValue() {
        // cheat.
        return oldValue;
    }

    public Dimension getDimension() {
        return Null.DIMENSION;
    }

    public String getLabel() {
        return "";
    }
}

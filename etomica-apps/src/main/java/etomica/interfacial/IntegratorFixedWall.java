/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.space.Space;

public class IntegratorFixedWall extends IntegratorVelocityVerlet {

    protected FixedWall fixedWall;
    
    public IntegratorFixedWall(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);
    }
    
    public void setFixedWall(FixedWall fixedWall) {
        this.fixedWall = fixedWall;
        getEventManager().addListener(fixedWall);
    }
    
    public void randomizeMomenta() {
        super.randomizeMomenta();
        if (fixedWall != null) fixedWall.integratorInitialized(null);
    }
}

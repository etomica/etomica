/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.kmc;

import etomica.space.Vector;
import etomica.simulation.Simulation;
import etomica.space.Space;

public class SimKMCworker extends Simulation{

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    public SimKMCworker(Space space) {
        super(space);
    }

    public static void main(String[] args){
        String fileName = args[0];
        final String APP_NAME = "SimKMCworker";

        final SimKMCMEAMadatom sim = new SimKMCMEAMadatom();
        Vector vect = sim.getSpace().makeVector();
        vect.setX(0, 9.8);
        vect.setX(1, -0.2);
        vect.setX(2, -0.2);
        
        sim.setMovableAtoms(100.0, vect);
        
        sim.setPotentialListAtoms();
        
        sim.initializeConfiguration("searchStart");
        sim.randomizePositions();
        
        sim.enableDimerSearch(fileName, 1000);
        sim.integratorDimer.setRotNum(0);
        
        sim.getController().actionPerformed();
    }
}

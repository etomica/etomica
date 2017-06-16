/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.kmc;

import etomica.space.Vector;
import etomica.simulation.Simulation;
import etomica.space.Space;

public class SimKMCmaster extends Simulation{

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    public SimKMCmaster(Space space) {
        super(space);
    }

    public static void main(String[] args){
        double temp = Double.parseDouble(args[0]);
        int steps = Integer.parseInt(args[1]);
        int totalSearch = Integer.parseInt(args[2]);
        final String APP_NAME = "SimKMCmaster";

        final SimKMCMEAMadatom sim = new SimKMCMEAMadatom();
        Vector vect = sim.getSpace().makeVector();
        vect.setX(0, 9.8);
        vect.setX(1, -0.2);
        vect.setX(2, -0.2);
               
        sim.setMovableAtoms(100.0, vect);
        
        sim.setPotentialListAtoms();
        
        sim.initializeConfiguration("initialStart");
        
        sim.integratorKMCCluster(temp, steps, totalSearch);
        
        //for sn energy: -3331480.584975273    Vib: 9.561284069712113E96
        sim.integratorKMCCluster.setInitialStateConditions(-3331480.5785103873, 8.148332378344422E277);
        sim.getController().actionPerformed();
    }
}

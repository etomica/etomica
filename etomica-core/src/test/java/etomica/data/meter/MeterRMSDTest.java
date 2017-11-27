/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by andrew on 6/15/17.
 */
public class MeterRMSDTest {
    Box box;
    MeterRMSD meter;
    double EPSILON = 1e-8;
    
    @Before
    public void setUp() throws Exception {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        box = new Box(space);
        sim.addBox(box);
        Vector L = space.makeVector();
        L.E(10);
        box.getBoundary().setBoxSize(L);
        Species species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        box.setNMolecules(species, 1);
        meter = new MeterRMSD(box, space);
    }
    
    @Test
    public void testGetDataAsScalar() throws Exception {
        Vector p = box.getLeafList().getAtom(0).getPosition();
        assertEquals(meter.getDataAsScalar(), 0, EPSILON);
        
        // move a full box length without apply PBC to atom
        for (double x = 0.5; x < 18.5; x++) {
            p.setX(0, x);
            assertEquals(meter.getDataAsScalar(), x, EPSILON);
        }
        
        // move atom back through the origin and more than a full box length
        // in the negative direction.  apply PBC each time.  PBC applied
        // first time causes atom to just 2 box lengths in a single step.
        for (double x = 18.5; x > -18.5; x--) {
            p.setX(0, x);
            Vector shift = box.getBoundary().centralImage(p);
            p.PE(shift);
            assertEquals(meter.getDataAsScalar(), Math.abs(x), EPSILON);
        }
    }
}
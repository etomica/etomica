/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorRigidIterative;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.molecule.OrientationCalcAtom;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotatingMolecule;

public class SingleDipole {

    public static SimulationGraphic makeDipole() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        Box box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
        sim.addBox(box);
        SpeciesSpheresRotatingMolecule species = new SpeciesSpheresRotatingMolecule(sim, space);
        sim.addSpecies(species);
        box.setNMolecules(species, 1);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.002;
        int maxIterations = 20;
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim, potentialMaster, timeInterval, 1, space);
//        integrator.printInterval = 10;
        integrator.setMaxIterations(maxIterations);
        integrator.setBox(box);
        OrientationCalcAtom calcer = new OrientationCalcAtom();
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(1);
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
        sim.getController().addAction(ai);

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

        ai.setSleepPeriod(10);
        SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1, space, sim.getController());
        
        return graphic;
    }
    
    public static void main(String[] args) {
        SimulationGraphic graphic = makeDipole();
        graphic.makeAndDisplayFrame();
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            SimulationGraphic graphic = makeDipole();

            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}

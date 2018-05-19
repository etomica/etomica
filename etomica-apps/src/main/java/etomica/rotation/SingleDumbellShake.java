/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerletShake;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.units.Kelvin;
import etomica.util.Constants;

public class SingleDumbellShake {

    public static SimulationGraphic makeSingleWater() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesSpheres species = new SpeciesSpheres(sim, space, 2);
        ((ConformationLinear) species.getConformation()).setAngle(0, 0);
        ((ConformationLinear) species.getConformation()).setAngle(1, 0.5 * Math.PI);
        ((ConformationLinear) species.getConformation()).setBondLength(2);
        sim.addSpecies(species);
        Box box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
        sim.addBox(box);
        box.setNMolecules(species, 1);
        box.setDensity(0.01 / 18.0 * Constants.AVOGADRO / 1E24);
//        new ConfigurationLattice(new LatticeCubicFcc(), space).initializeCoordinates(box);
        IMolecule molecule = box.getMoleculeList(species).get(0);
        species.getConformation().initializePositions(molecule.getChildList());
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeStep = 2 * Math.PI / 100;
        int maxIterations = 20;
        IntegratorVelocityVerletShake integrator = new IntegratorVelocityVerletShake(sim, potentialMaster, box);
        integrator.setTimeStep(timeStep);
//        integrator.printInterval = 10;
        integrator.setMaxIterations(maxIterations);
        integrator.setBondConstraints(species, new int[][]{{0, 1}}, new double[]{2});
        integrator.setIsothermal(false);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
//        integrator.setThermostatInterval(100);
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
//        System.out.println("using rigid with dt="+dt);
        sim.getController().addAction(ai);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

        if (false) {
            ai.setMaxSteps(100);
            sim.getController().actionPerformed();
            return null;
        }
        ai.setSleepPeriod(10);
        SimulationGraphic graphic = new SimulationGraphic(sim, "SHAKE", 1);
        return graphic;
    }
    
    public static void main(String[] args) {
        SimulationGraphic graphic = makeSingleWater();
        if (graphic != null) {
            graphic.makeAndDisplayFrame();
        }
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            SimulationGraphic graphic = makeSingleWater();

            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}

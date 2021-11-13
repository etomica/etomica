/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOrientedKinetic;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorRigidIterative;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.OrientationCalcAtomLinear;
import etomica.potential.P2LJDipoleAtomic;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Pixel;

public class DipoleBox extends Simulation {

    public final IntegratorRigidIterative integrator;
    public final Box box;

    public DipoleBox(Space space, int nAtoms, double dt) {
        super(space);
        SpeciesGeneral species = SpeciesSpheresRotating
                .create(space, new ElementSimple("A"), true, true);
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularPeriodic(getSpace(), 10));
        box.setNMolecules(species, nAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.5);
        inflater.actionPerformed();
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < nAtoms; i++) {
            IAtomOrientedKinetic atom = (IAtomOrientedKinetic) atoms.get(i);
            IOrientation orientation = atom.getOrientation();
            for (int j = 0; j < 20; j++) {
                orientation.randomRotation(getRandom(), 1);
            }
        }
        NeighborManagerSimple neighborManager = new NeighborManagerSimple(box);
        PotentialComputePairGeneral potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);
        int maxIterations = 20;
        integrator = new IntegratorRigidIterative(getSpeciesManager(), random, potentialMaster, dt, 1, box);
        integrator.setTemperature(1);
        integrator.setIsothermal(false);
        integrator.printInterval = 100;
        integrator.setMaxIterations(maxIterations);
        OrientationCalcAtomLinear calcer = new OrientationCalcAtomLinear();
        integrator.setOrientationCalc(species, calcer);
        this.getController().addActivity(new ActivityIntegrate(integrator));

        P2LJDipoleAtomic p2 = new P2LJDipoleAtomic(space, 1.0, 1.0, 2.0);
        p2.setTruncationRadius(2.5);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2);

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        int nAtoms = 864;
        double dt = 0.01;
        if (args.length == 0) {
            DipoleBox sim = new DipoleBox(space, nAtoms, dt);
            sim.getController().setSleepPeriod(10);
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
            graphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(30));
            graphic.makeAndDisplayFrame();
        }
        else {
            dt = Double.parseDouble(args[0]);
            nAtoms = Integer.parseInt(args[1]);
            DipoleBox sim = new DipoleBox(space, nAtoms, dt);
            sim.integrator.printInterval = 100;
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
        }
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            Space space = Space3D.getInstance();
            DipoleBox sim = new DipoleBox(space, 864, 0.01);
            sim.getController().setSleepPeriod(10);
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
            graphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(30));

            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}

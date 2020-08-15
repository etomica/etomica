/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorRigidMatrixIterative;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.*;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresRotatingMolecule;
import etomica.units.Pixel;

public class DipoleBoxMatrix extends Simulation {

    public final IntegratorRigidMatrixIterative integrator;
    public final Box box;
    
    public DipoleBoxMatrix(Space space, int nAtoms, double dt) {
        super(space);
        SpeciesSpheresRotatingMolecule species = new SpeciesSpheresRotatingMolecule(this, space, Vector.of(new double[]{0.025, 0.025, 0.025}));
        species.setIsDynamic(true);
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularPeriodic(getSpace(), 10));
        box.setNMolecules(species, nAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.5);
        inflater.actionPerformed();
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        IMoleculeList molecules = box.getMoleculeList();
        for (int i = 0; i < nAtoms; i++) {
            IMolecule molecule = molecules.get(i);
            ((IMoleculePositioned) molecule).getPosition().E(molecule.getChildList().get(0).getPosition());
            IOrientationFull3D orientation = (IOrientationFull3D) ((IMoleculeOriented) molecule).getOrientation();
            for (int j = 0; j < 20; j++) {
                orientation.randomRotation(getRandom(), 1);
            }
        }
        PotentialMaster potentialMaster = new PotentialMaster();
        int maxIterations = 20;
        integrator = new IntegratorRigidMatrixIterative(this, potentialMaster, dt, 1, box);
        integrator.setTemperature(1);
        integrator.setIsothermal(false);
        integrator.printInterval = 1;
        integrator.setMaxIterations(maxIterations);
        OrientationCalcAtom calcer = new OrientationCalcAtom();
        integrator.setOrientationCalc(species, calcer);
        this.getController().addActivity(new ActivityIntegrate(integrator));

        P2LJDipole p2 = new P2LJDipole(space, 1.0, 1.0, 2.0);
        p2.setTruncationRadius(2.5);
        potentialMaster.addPotential(p2, new ISpecies[]{species, species});

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        int nAtoms = 864;
        double dt = 0.01;
        if (args.length == 0) {
            DipoleBoxMatrix sim = new DipoleBoxMatrix(space, nAtoms, dt);
            sim.getController().setSleepPeriod(10);
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
            graphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(30));
            graphic.makeAndDisplayFrame();
        }
        else {
            dt = Double.parseDouble(args[0]);
            nAtoms = Integer.parseInt(args[1]);
            DipoleBoxMatrix sim = new DipoleBoxMatrix(space, nAtoms, dt);
            sim.integrator.printInterval = 100;
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), Long.MAX_VALUE);
        }
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            Space space = Space3D.getInstance();
            DipoleBoxMatrix sim = new DipoleBoxMatrix(space, 864, 0.01);
            sim.getController().setSleepPeriod(10);
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
            graphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(30));

            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}

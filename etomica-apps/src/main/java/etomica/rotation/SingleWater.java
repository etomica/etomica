/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.BoxImposePbc;
import etomica.action.IAction;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorRigidIterative;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.SpeciesWater3POriented;
import etomica.molecule.MoleculeOriented;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.random.RandomNumberGenerator;

import java.awt.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SingleWater {

    public static SimulationGraphic makeSingleWater() {
        final Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        ((RandomNumberGenerator) sim.getRandom()).getWrappedRandom().setSeed(2000);
        SpeciesWater3POriented species = new SpeciesWater3POriented(sim.getSpace(), true);
        sim.addSpecies(species);
        final Box box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
        sim.addBox(box);
        box.setNMolecules(species, 1);
        box.setDensity(0.01 / 18.0 * Constants.AVOGADRO / 1E24);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.00016;
        int maxIterations = 40;
        final IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim, potentialMaster, timeInterval, 1, box);
        integrator.printInterval = 100;
        integrator.setMaxIterations(maxIterations);
        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
//        integrator.setOrientAtom((IAtomPositioned)((IMolecule)box.getMoleculeList(speciesOrient).getAtom(0)).getChildList().getAtom(0));
        integrator.setOrientationCalc(species, calcer);
//        integrator.setIsothermal(true);
        integrator.setTemperature(Kelvin.UNIT.toSim(148.5));
//        integrator.setThermostatInterval(100);
        sim.getController().addActivity(new ActivityIntegrate(integrator));
//        System.out.println("using rigid with dt="+dt);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

        if (false) {
            final boolean isWriting = false;
            final FileWriter fileWriter;
            FileReader fileReader;
            final BufferedReader bufReader;
            if (isWriting) {
                try {
                    fileWriter = new FileWriter("rigidA008.out");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                bufReader = null;
            } else {
                fileWriter = null;
                String infile = "rigidA008.out";
                try {
                    fileReader = new FileReader(infile);
                } catch (IOException e) {
                    throw new RuntimeException("Cannot open " + infile + ", caught IOException: " + e.getMessage());
                }
                bufReader = new BufferedReader(fileReader);
            }

            IAction writeA = new IAction() {
                public void actionPerformed() {
                    MoleculeOriented molecule = (MoleculeOriented) box.getMoleculeList().get(0);
                    A.setOrientation((IOrientationFull3D) molecule.getOrientation());
                    try {
                        if (isWriting) {
                            fileWriter.write(integrator.getCurrentTime() + " ");
                            for (int i = 0; i < 3; i++) {
                                for (int j = 0; j < 3; j++) {
                                    fileWriter.write(A.component(i, j) + " ");
                                }
                            }
                            fileWriter.write("\n");
                        } else {
                            String line = bufReader.readLine();
                            String[] componentsStr = line.split(" +");
                            int k = 1;
                            for (int i = 0; i < 3; i++) {
                                for (int j = 0; j < 3; j++) {
                                    Aex.setComponent(i, j, Double.parseDouble(componentsStr[k]));
                                    k++;
                                }
                            }
                            A.ME(Aex);
                            Aex.E(A);
                            Aex.transpose();
                            A.TE(Aex);
                            double err = A.trace() / 6.0;
                            sum += err;
                            n++;
                            System.out.println(integrator.getCurrentTime() + " " + Math.sqrt(sum / n));
                        }
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }

                double sum = 0;
                int n = 0;
                RotationTensor3D A = (RotationTensor3D) space.makeRotationTensor();
                RotationTensor3D Aex = (RotationTensor3D) space.makeRotationTensor();
            };
            IntegratorListenerAction writeAListener = new IntegratorListenerAction(writeA);
            writeAListener.setInterval(100);
            integrator.getEventManager().addListener(writeAListener);
            sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, Long.MAX_VALUE));
        } else {
//          ai.setSleepPeriod(10);
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getOxygenType(), Color.RED);
            return graphic;
        }
        return null;
    }
    
    public static void main(String[] args) {
        SimulationGraphic graphic = makeSingleWater();
        graphic.makeAndDisplayFrame();
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            SimulationGraphic graphic = makeSingleWater();

            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }
}

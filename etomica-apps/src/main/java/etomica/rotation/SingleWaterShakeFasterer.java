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
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.integrator.ShakeListener;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.water.ConformationWater3P;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.SpeciesWater3P;
import etomica.molecule.IMolecule;
import etomica.potential.BondingInfo;
import etomica.potential.PotentialMasterFasterer;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesAgentManager;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SingleWaterShakeFasterer {

    public static SimulationGraphic makeSingleWater() {
        final Space space = Space3D.getInstance();
        final Simulation sim = new Simulation(space);
        SpeciesGeneral species = SpeciesWater3P.create(true);
        sim.addSpecies(species);
        final Box box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
        sim.addBox(box);
        box.setNMolecules(species, 1);
        box.setDensity(0.01 / 18.0 * Constants.AVOGADRO / 1E24);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        double timeStep = 0.000166;
        int maxIterations = 200;
        PotentialMasterFasterer potentialMaster = new PotentialMasterFasterer(sim.getSpeciesManager(), box, BondingInfo.noBonding());
        final IntegratorVelocityVerletFasterer integrator = new IntegratorVelocityVerletFasterer(potentialMaster, sim.getRandom(), timeStep, Kelvin.UNIT.toSim(271.2654804973), box);
        integrator.setTimeStep(timeStep);
        double lOH = ConformationWater3P.bondLengthOH;
        double lHH = Math.sqrt(2 * lOH * lOH * (1 - Math.cos(ConformationWater3P.angleHOH)));
        SpeciesAgentManager<ShakeListener.BondConstraints> shakeAgents = new SpeciesAgentManager<>(new SpeciesAgentManager.AgentSource<ShakeListener.BondConstraints>() {
            @Override
            public ShakeListener.BondConstraints makeAgent(ISpecies type) {
                return new ShakeListener.BondConstraints(new int[][]{{0, 2}, {1, 2}, {0, 1}}, new double[]{lOH, lOH, lHH});
            }

            @Override
            public void releaseAgent(ShakeListener.BondConstraints agent, ISpecies type) {

            }
        }, sim.getSpeciesManager());
        ShakeListener shake = new ShakeListener(sim.getSpeciesManager(), shakeAgents, integrator);
        shake.setMaxIterations(maxIterations);
        integrator.getEventManager().addListener(shake);
        integrator.setIsothermal(false);
//        integrator.setThermostatInterval(100);
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
                    fileWriter = new FileWriter("shakeA008.out");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                bufReader = null;
            } else {
                fileWriter = null;
                String infile = "shakeA008.out";
                try {
                    fileReader = new FileReader(infile);
                } catch (IOException e) {
                    throw new RuntimeException("Cannot open " + infile + ", caught IOException: " + e.getMessage());
                }
                bufReader = new BufferedReader(fileReader);
            }

            final OrientationCalcWater3P calcer = new OrientationCalcWater3P(space);
            final IOrientationFull3D orientation = (IOrientationFull3D) calcer.makeOrientation(space);
            IAction writeA = new IAction() {
                public void actionPerformed() {
                    IMolecule molecule = box.getMoleculeList().get(0);
                    calcer.calcOrientation(molecule, orientation);
                    A.setOrientation(orientation);
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
            sim.getController().setSleepPeriod(1);
            sim.getController().addActivity(new ActivityIntegrate(integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim, "SHAKE", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("H"), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("O"), Color.RED);
            return graphic;
        }
        return null;
    }

    public static void main(String[] args) {
        SimulationGraphic graphic = makeSingleWater();
        if (graphic != null) {
            graphic.makeAndDisplayFrame();
        }
    }

}

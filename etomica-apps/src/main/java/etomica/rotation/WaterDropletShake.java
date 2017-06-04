/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerletShake;
import etomica.models.water.ConformationWater3P;
import etomica.models.water.P2WaterSPC;
import etomica.models.water.SpeciesWater3P;
import etomica.potential.P2Electrostatic;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.awt.*;

public class WaterDropletShake {

    public static SimulationGraphic makeWaterDroplet() {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        Box box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        sim.addBox(box);
        SpeciesWater3P species = new SpeciesWater3P(sim.getSpace(), true);
        sim.addSpecies(species);
        box.setNMolecules(species, 108);
        box.setDensity(1/18.0*Constants.AVOGADRO/1E24);
        ConfigurationWater108 configFile = new ConfigurationWater108();
        configFile.initializeCoordinates(box);
        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.002;
        int maxIterations = 100;
        IntegratorVelocityVerletShake integrator = new IntegratorVelocityVerletShake(sim, potentialMaster, space);
        double lOH = ConformationWater3P.bondLengthOH;
        double lHH = Math.sqrt(2*lOH*lOH*(1-Math.cos(ConformationWater3P.angleHOH)));
        integrator.setBondConstraints(species, new int[][]{{0,2},{1,2},{0,1}}, new double[]{lOH, lOH, lHH});
        integrator.setTimeStep(timeInterval);
        integrator.setMaxIterations(maxIterations);
        integrator.setBox(box);
//        integrator.setOrientAtom((IAtomPositioned)((IMolecule)box.getMoleculeList(speciesOrient).getAtom(0)).getChildList().getAtom(0));
        integrator.setIsothermal(false);
        integrator.setTemperature(Kelvin.UNIT.toSim(273));
        try {
            integrator.reset();
        }
        catch (ConfigurationOverlapException e){}
//        MeterTemperature meterTemperature = new MeterTemperature(box, 3);
//        System.out.println("T="+Kelvin.UNIT.fromSim(meterTemperature.getDataAsScalar()));
//        integrator.setThermostatInterval(100);
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
//        System.out.println("using rigid with dt="+dt);
        sim.getController().addAction(ai);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

        
        double chargeOxygen = Electron.UNIT.toSim(-0.82);
        double chargeHydrogen = Electron.UNIT.toSim(0.41);

        AtomType oType = species.getOxygenType();
        AtomType hType = species.getHydrogenType();
        double epsOxygen = new P2WaterSPC(space).getEpsilon();
        double sigOxygen = new P2WaterSPC(space).getSigma();
        PotentialGroup pGroup = potentialMaster.makePotentialGroup(2);
        P2LennardJones potentialLJOO = new P2LennardJones(space, sigOxygen, epsOxygen);
        pGroup.addPotential(potentialLJOO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{oType, oType}));

        P2Electrostatic potentialQHH = new P2Electrostatic(space);
        potentialQHH.setCharge1(chargeHydrogen);
        potentialQHH.setCharge2(chargeHydrogen);
        pGroup.addPotential(potentialQHH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{hType, hType}));

        P2Electrostatic potentialQOO = new P2Electrostatic(space);
        potentialQOO.setCharge1(chargeOxygen);
        potentialQOO.setCharge2(chargeOxygen);
        pGroup.addPotential(potentialQOO, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{oType, oType}));
        
        P2Electrostatic potentialQOH = new P2Electrostatic(space);
        potentialQOH.setCharge1(chargeOxygen);
        potentialQOH.setCharge2(chargeHydrogen);
        pGroup.addPotential(potentialQOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{oType, hType}));
        pGroup.addPotential(potentialQOH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{hType, oType}));

        potentialMaster.addPotential(pGroup, new ISpecies[]{species,species});

        if (false) {
            ai.setSleepPeriod(2);
            SimulationGraphic graphic = new SimulationGraphic(sim, "SHAKE", 1, space, sim.getController());
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getOxygenType(), Color.RED);
            return graphic;
        }
        sim.getController().actionPerformed();
        return null;
    }

    public static void main(String[] args) {
        SimulationGraphic graphic = makeWaterDroplet();
        if (graphic != null) {
            graphic.makeAndDisplayFrame();
        }
    }
    
    public static class Applet extends javax.swing.JApplet {

        private static final long serialVersionUID = 1L;

        public void init() {
            SimulationGraphic graphic = makeWaterDroplet();

            getContentPane().add(graphic.getPanel());
        }
    }
}

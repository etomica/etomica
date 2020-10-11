/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.action.BoxImposePbc;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionAll;
import etomica.potential.EwaldSummation;
import etomica.potential.EwaldSummation.MyCharge;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Bar;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import java.awt.*;
import java.util.ArrayList;

public class TestEwaldTIP4PWater extends Simulation {


    private final static String APP_NAME = "Test Ewald Sum TIP4P Water";
    private static final int PIXEL_SIZE = 15;
    protected final PotentialMaster potentialMaster;
    protected final Box box;
    protected final SpeciesGeneral species;
    protected final IntegratorMC integrator;

    TestEwaldTIP4PWater(Space space) {
        super(space);

        LatticeCubicFcc lattice = new LatticeCubicFcc(space);
        ConfigurationLattice configuration = new ConfigurationLattice(lattice, space);

        ConformationWaterTIP4P config = new ConformationWaterTIP4P(space);
        species = SpeciesWater4P.create(config);
        addSpecies(species);

        potentialMaster = new PotentialMaster();

        box = this.makeBox();
        box.getBoundary().setBoxSize(Vector.of(new double[]{25, 25, 25}));
        box.setNMolecules(species, 125);

        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));

        MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(this, potentialMaster, space);
        MCMoveRotateMolecule3D mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random, space);
        MCMoveVolume mcMoveVolume = new MCMoveVolume(potentialMaster, random, space, Bar.UNIT.toSim(1.0132501));

        ((MCMoveStepTracker) mcMoveVolume.getTracker()).setNoisyAdjustment(true);


        integrator.getMoveManager().addMCMove(mcMoveMolecule);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        integrator.getMoveManager().addMCMove(mcMoveVolume);


        getController().addActivity(new ActivityIntegrate(integrator), 6000);


        //Potential
        P2LennardJones potentialLJ = new P2LennardJones(space, 3.154, Kelvin.UNIT.toSim(78.02));
        potentialMaster.addPotential(potentialLJ, new AtomType[]{species.getTypeByName("O"), species.getTypeByName("O")});

        CriterionAll criterionAll = new CriterionAll();

        //Ewald Summation
        ChargeAgentSourceTIP4PWater agentSource = new ChargeAgentSourceTIP4PWater();
        AtomLeafAgentManager<MyCharge> atomAgentManager = new AtomLeafAgentManager<MyCharge>(agentSource, box);
        EwaldSummation ewaldSummation = new EwaldSummation(box, atomAgentManager, space, 4, 9);
//		ewaldSummation.setCriterion(criterionAll);
//		ewaldSummation.setBondedIterator(new ApiIntragroup());
        potentialMaster.addPotential(ewaldSummation, new AtomType[0]);
        ////////


        BoxImposePbc imposePBC = new BoxImposePbc(box, space);

        configuration.initializeCoordinates(box);

        integrator.getEventManager().addListener(new IntegratorListenerAction(imposePBC));

    }

    public static void main (String[] args){

        Space sp = Space3D.getInstance();
		TestEwaldTIP4PWater sim = new TestEwaldTIP4PWater(sp);
		SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 1);
		Pixel pixel = new Pixel(10);
		simGraphic.getDisplayBox(sim.box).setPixelUnit(pixel);
		ArrayList dataStreamPumps = simGraphic.getController().getDataStreamPumps();

        /////////////////////////////////////////////////////////////
		MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
		DisplayTextBox PEbox = new DisplayTextBox();
		DataPump PEpump = new DataPump(meterPE, PEbox);
		dataStreamPumps.add(PEpump);

        IntegratorListenerAction pumpListener = new IntegratorListenerAction(PEpump);
		pumpListener.setInterval(1);
	    sim.integrator.getEventManager().addListener(pumpListener);

        simGraphic.add(PEbox);
	    //////////////////////////////////////////////////////////


        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getTypeByName("H"), Color.WHITE);
        colorScheme.setColor(sim.species.getTypeByName("O"), Color.RED);
        ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getTypeByName("M"), 0);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        simGraphic.getDisplayBox(sim.box).repaint();

	}
}

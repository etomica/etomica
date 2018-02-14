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
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.CriterionAll;
import etomica.potential.EwaldSummation;
import etomica.potential.EwaldSummation.MyCharge;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Bar;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import java.util.ArrayList;

public class TestEwaldTIP4PWater extends Simulation {


    private final static String APP_NAME = "Test Ewald Sum TIP4P Water";
    private static final int PIXEL_SIZE = 15;
    protected final PotentialMaster potentialMaster;
    protected final Box box;
    protected final SpeciesWater4P species;
    protected final IntegratorMC integrator;
    protected final BoundaryRectangularPeriodic boundary;

    TestEwaldTIP4PWater(Space space){
		super(space);
		potentialMaster = new PotentialMaster();

		LatticeCubicFcc lattice = new LatticeCubicFcc(space);
		ConfigurationLattice configuration = new ConfigurationLattice(lattice, space);

		ConformationWaterTIP4P config = new ConformationWaterTIP4P(space);
		species = new SpeciesWater4P(space);
		species.setConformation(config);
		addSpecies(species);

		box = new Box(space);
		addBox(box);
		box.getBoundary().setBoxSize(space.makeVector(new double[] {25, 25, 25}));
		box.setNMolecules(species, 125);

		integrator = new IntegratorMC(this, potentialMaster, box);
		integrator.setTemperature(Kelvin.UNIT.toSim(298));

		MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(this, potentialMaster, space);
		MCMoveRotateMolecule3D mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random, space);
		MCMoveVolume mcMoveVolume = new MCMoveVolume(potentialMaster, random, space, Bar.UNIT.toSim(1.0132501));

		((MCMoveStepTracker)mcMoveVolume.getTracker()).setNoisyAdjustment(true);


        integrator.getMoveManager().addMCMove(mcMoveMolecule);
		integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
		integrator.getMoveManager().addMCMove(mcMoveVolume);


        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(6000);
        getController().addAction(activityIntegrate);



        //Potential
		P2LennardJones potentialLJ = new P2LennardJones(space, 3.154,Kelvin.UNIT.toSim(78.02));
        potentialMaster.addPotential(potentialLJ, new AtomType[]{species.getOxygenType(), species.getOxygenType()});

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

        boundary = new BoundaryRectangularPeriodic(space, 20);
        boundary.setBoxSize(space.makeVector(new double[] {20, 20, 20}));
        box.setBoundary(boundary);

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
		MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
		meterPE.setBox(sim.box);
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
        colorScheme.setColor(sim.species.getHydrogenType(), java.awt.Color.WHITE);
        colorScheme.setColor(sim.species.getOxygenType(), java.awt.Color.RED);
        ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getMType(), 0);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        simGraphic.getDisplayBox(sim.box).repaint();

	}
}

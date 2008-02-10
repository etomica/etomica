package etomica.models.water;

import java.util.ArrayList;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.ApiIntragroup;
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
import etomica.nbr.CriterionAll;
import etomica.potential.EwaldSummation;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.units.Bar;
import etomica.units.Kelvin;
import etomica.units.Pixel;

public class TestEwaldTIP4PWater extends Simulation {
	
	
	TestEwaldTIP4PWater(){
		super(Space3D.getInstance(), true);
		potentialMaster = new PotentialMaster(space);
		
		LatticeCubicFcc lattice = new LatticeCubicFcc();
		ConfigurationLattice configuration = new ConfigurationLattice(lattice);
		
		ConformationWaterTIP4P config = new ConformationWaterTIP4P(space);
		species = new SpeciesWater4P(space);
		species.getMoleculeType().setConformation(config);
		getSpeciesManager().addSpecies(species);
		
		integrator = new IntegratorMC(this, potentialMaster);
		integrator.setTemperature(Kelvin.UNIT.toSim(298));
		
		MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(this, potentialMaster);
		MCMoveRotateMolecule3D mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
		MCMoveVolume mcMoveVolume = new MCMoveVolume(potentialMaster, random, Bar.UNIT.toSim(1.0132501));
		
		((MCMoveStepTracker)mcMoveVolume.getTracker()).setNoisyAdjustment(true);
		
		
		integrator.getMoveManager().addMCMove(mcMoveMolecule);
		integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
		integrator.getMoveManager().addMCMove(mcMoveVolume);
		
		
		
		
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(6000);
        getController().addAction(activityIntegrate);
        
		box = new Box(this);
		addBox(box);
		box.setDimensions(Space.makeVector(new double[] {25, 25, 25}));
		box.setNMolecules(species, 125);
		
		
		//Potential
		P2LennardJones potentialLJ = new P2LennardJones(space, 3.154,Kelvin.UNIT.toSim(78.02));
		potentialMaster.addPotential(potentialLJ, new AtomType[]{species.getOxygenType(), species.getOxygenType()} );
        
		CriterionAll criterionAll = new CriterionAll();
		
		//Ewald Summation
		ChargeAgentSourceTIP4PWater agentSource = new ChargeAgentSourceTIP4PWater();
		AtomAgentManager atomAgentManager = new AtomAgentManager(agentSource, box);
		EwaldSummation ewaldSummation = new EwaldSummation(box, atomAgentManager, 0);
		ewaldSummation.setCriterion(criterionAll);
		ewaldSummation.setBondedIterator(new ApiIntragroup());
		potentialMaster.addPotential(ewaldSummation, new Species[0]);
		////////
		
		
		
		BoxImposePbc imposePBC = new BoxImposePbc(box);
		
        boundary = new BoundaryRectangularPeriodic(space, getRandom(), 20);
        boundary.setDimensions(Space.makeVector(new double[] {20, 20, 20}));
        box.setBoundary(boundary);
        
        configuration.initializeCoordinates(box);
        
        integrator.setBox(box);
        integrator.addIntervalAction(imposePBC);
		
	}
	
	public static void main (String[] args){
		
		TestEwaldTIP4PWater sim = new TestEwaldTIP4PWater();
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
		
	    sim.integrator.addIntervalAction(PEpump);
	    sim.integrator.setActionInterval(PEpump, 1);
		
	    simGraphic.add(PEbox);
	    //////////////////////////////////////////////////////////
	    
	    
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        AtomTypeMolecule atomType = sim.species.getMoleculeType();
        colorScheme.setColor(atomType.getChildTypes()[0], java.awt.Color.white);
        colorScheme.setColor(atomType.getChildTypes()[1], java.awt.Color.blue);
        ((AtomTypeSphere)atomType.getChildTypes()[2]).setDiameter(0);
        
        simGraphic.makeAndDisplayFrame(APP_NAME);

        simGraphic.getDisplayBox(sim.box).repaint();

	}
	
	protected final PotentialMaster potentialMaster;
	protected final Box box;
	protected final SpeciesWater4P species;
	protected final IntegratorMC integrator;
	protected final BoundaryRectangularPeriodic boundary;
	private final static String APP_NAME = "Test Ewald Sum TIP4P Water";
	private static final int PIXEL_SIZE = 15;
	private static final long serialVersionUID = 1L;
}

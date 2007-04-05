package etomica.yukawa;

import java.awt.Color;

import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

/**
 * A Yukawa MD Simulation
 * 
 * adapted from HSMD3D simulation prototype.  employs neighbor listing.
 * 
 * @author msellers
 *
 */

public class TestYukawaMD3D extends Simulation{
	public final Phase phase;
	public final IntegratorVelocityVerlet integrator;
	public final SpeciesSpheresMono species;
	public final P2Yukawa potential;
	
	public TestYukawaMD3D(){
		this(new Default());
	}
	
	public TestYukawaMD3D(Default defaults){
		this(Space3D.getInstance(), defaults);
	}
	
	private TestYukawaMD3D(Space space, Default defaults){
		super(space, true, new PotentialMasterList(space, 1.6), Default.BIT_LENGTH, defaults);
		
		int numAtoms = 256;
		double neighborRangeFac = 1.6;
		
		defaults.makeLJDefaults();
		defaults.atomSize = 1.0;
		defaults.boxSize = 14.4573*Math.pow((numAtoms/2020.0),1.0/3.0);
		
		((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*defaults.atomSize);
		
		integrator = new IntegratorVelocityVerlet(this);
		integrator.setIsothermal(false);
		integrator.setTimeStep(0.01);
		this.register(integrator);
		
		
		ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
		activityIntegrate.setDoSleep(true);
		activityIntegrate.setSleepPeriod(1);
		getController().addAction(activityIntegrate);
		
		species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
		phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);
        NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager(phase);
        integrator.addListener(nbrManager);
		potential = new P2Yukawa(this);
		
		double truncationRadius = 2.5*potential.getKappa();
		if(truncationRadius > 0.5*phase.getBoundary().getDimensions().x(0)){
			throw new RuntimeException("Truncaiton radius too large.  Max allowed is "+0.5*phase.getBoundary().getDimensions().x(0));
		}
		P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, truncationRadius);
		((PotentialMasterList)potentialMaster).setCellRange(3);
		((PotentialMasterList)potentialMaster).setRange(potentialTruncated.getRange()*1.3);
		potentialMaster.addPotential(potentialTruncated, new Species[] {species, species});
		
		new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(phase);
		integrator.setPhase(phase);
	}
	
	public static void main(String[] args){
		Default defaults = new Default();
		defaults.doSleep = false;
		defaults.ignoreOverlap = true;
		
		TestYukawaMD3D sim = new TestYukawaMD3D(defaults);
		SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpeciesAgent(sim.phase.getAgent(sim.species));
		simGraphic.add(nSelector);
		simGraphic.makeAndDisplayFrame();
		ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
		colorScheme.setColor(sim.species.getMoleculeType(), Color.red);
	}
}

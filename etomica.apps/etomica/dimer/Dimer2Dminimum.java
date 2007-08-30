package etomica.dimer;

import etomica.action.BoxImposePbc;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation using Henkelman's Dimer method to find PE minimum from saddle point.
 * 
 * @author msellers
 *
 */

public class Dimer2Dminimum extends Simulation {

	private static final long serialVersionUID = 1L;
	public ActivityIntegrate activityIntegrateDimerRT, activityIntegrateDimerMEP;
	public IntegratorDimerMEP integratorDimerMEP;
	public SpeciesSpheresMono species;
	public Box box;
	public P1LepsHarmonic potential;
	public PotentialMaster potentialMaster;
	public MeterPotentialEnergy energy;
	public AccumulatorAverage avgEnergy;
	public DataPump pump;
	
	public static void main(String[] args){
		final String APP_NAME = "Dimer2Dminumim";
		final Dimer2Dminimum sim = new Dimer2Dminimum();

		//sim.activityIntegrateDimerMEP.setMaxSteps(500);

		// Write initial configuration
		WriteConfiguration writeConfig = new WriteConfiguration();
		writeConfig.setConfName("dimer-config-min");
		writeConfig.setBox(sim.box);
		writeConfig.actionPerformed();

		MeterPotentialEnergy energy = new MeterPotentialEnergy(sim.potentialMaster);
		
		energy.setBox(sim.box);

		DataLogger elog = new DataLogger();
		elog.setFileName("dimer-energy-min");
		elog.setAppending(true);
		elog.setCloseFileEachTime(true);
		elog.setDataSink(new DataTableWriter());
		sim.getController().getEventManager().addListener(elog);
		
		DataPump pump = new DataPump(energy, elog);
		
		sim.integratorDimerMEP.addIntervalAction(pump);
		sim.integratorDimerMEP.setActionInterval(pump, 1);
		
		sim.getController().actionPerformed();
		
		System.out.println("X Location: "+((IAtomPositioned)sim.box.getLeafList().getAtom(0)).getPosition().x(0));
		System.out.println("Y Location: "+((IAtomPositioned)sim.box.getLeafList().getAtom(0)).getPosition().x(1));
	}
		
	public Dimer2Dminimum() {
		
		super(Space2D.getInstance());
		potentialMaster = new PotentialMaster(space);
		
		species = new SpeciesSpheresMono(this);
		getSpeciesManager().addSpecies(species);
		
		box = new Box(this);
		addBox(box);
		box.setNMolecules(species, 1);
		
		potential = new P1LepsHarmonic(space);
		potentialMaster.addPotential(potential, new Species[]{species});

		BoxImposePbc imposepbc = new BoxImposePbc();
		imposepbc.setBox(box);
		
		// Create saddle point array
		IVector [] saddle = new IVector[1];
		saddle[0] = box.getSpace().makeVector();
		saddle[0].setX(0, 1.9348841569069544);
		saddle[0].setX(1, -1.2945894438743524);
		
		integratorDimerMEP = new IntegratorDimerMEP(this, potentialMaster, saddle, -10E-4, -0.05, new Species[]{species});
		
		activityIntegrateDimerMEP = new ActivityIntegrate(integratorDimerMEP);
		getController().addAction(activityIntegrateDimerMEP);
	}
}

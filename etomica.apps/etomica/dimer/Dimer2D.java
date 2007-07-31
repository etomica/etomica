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
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation using Henkelman's Dimer method.
 * 
 * @author msellers
 *
 */

public class Dimer2D extends Simulation {

	private static final long serialVersionUID = 1L;
	public ActivityIntegrate activityIntegrateDimer;
	public IntegratorDimerRT integratorDimer;
	public SpeciesSpheresMono species;
	public Box box;
	public P1LepsHarmonic potential;
	public PotentialMaster potentialMaster;
	public MeterPotentialEnergy energy;
	public AccumulatorAverage avgEnergy;
	public DataPump pump;
	
	public static void main(String[] args){
		final String APP_NAME = "Dimer2D";
		final Dimer2D sim = new Dimer2D();

		sim.activityIntegrateDimer.setMaxSteps(1000);

		// Write initial configuration
		WriteConfiguration writeConfig = new WriteConfiguration();
		writeConfig.setConfName("dimer-config");
		writeConfig.setBox(sim.box);
		writeConfig.actionPerformed();
		
		MeterPotentialEnergy energy = new MeterPotentialEnergy(sim.potentialMaster);
		
		energy.setBox(sim.box);

		DataLogger elog = new DataLogger();
		elog.setFileName("00_sim-energy");
		elog.setAppending(true);
		elog.setCloseFileEachTime(true);
		elog.setDataSink(new DataTableWriter());
		sim.getController().getEventManager().addListener(elog);
		
		DataPump pump = new DataPump(energy, elog);
		
		sim.integratorDimer.addIntervalAction(pump);
		sim.integratorDimer.setActionInterval(pump, 1);
		
		sim.getController().actionPerformed();
		
		
		System.out.println("---- Saddle Point Data ----");
		System.out.println("Steps: "+sim.integratorDimer.counter);
		System.out.println("Dimer Energy: "+energy.getDataAsScalar());
		System.out.println("Dimer Force Array Magnitude: "+sim.integratorDimer.saddleT);
		System.out.println("Dimer Rotational Force Magnitude: "+sim.integratorDimer.Frot);
		
		double Fmag = 0;
		for(int i=0;i<sim.integratorDimer.Feff.length;i++){
			Fmag += sim.integratorDimer.Feff[i].squared();
		}
		Fmag = Math.sqrt(Fmag);
		
		System.out.println("Dimer Translational Force Magnitude: "+Fmag);
		
		System.out.println("X Location: "+((IAtomPositioned)sim.box.getLeafList().getAtom(0)).getPosition().x(0));
		System.out.println("Y Location: "+((IAtomPositioned)sim.box.getLeafList().getAtom(0)).getPosition().x(1));
	}
		
	public Dimer2D() {
		
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
		
		double startx = 0.8;
		double starty = 0.8;
		
		((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition().setX(0, startx);
		((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition().setX(1, starty);
		
		
		
		integratorDimer = new IntegratorDimerRT(this, potentialMaster);
		
		activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
		getController().addAction(activityIntegrateDimer);
		
	}
}

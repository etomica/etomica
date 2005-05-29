package etomica.virial.simulations;

import etomica.*;
import etomica.space3d.Boundary;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.*;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorList;
import etomica.data.DataSourceCountSteps;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.MCMove;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.virial.*;
import etomica.virial.cluster.*;
import etomica.virial.simulations.*;
import etomica.models.water.*;
import etomica.modifier.ModifierBoolean;

/**
 * @author kofke
 *
 * Calculation of virial coefficients for water models.
 */
public class SimulationVirialWater extends SimulationGraphic {

	/**
	 * Default constructor, using a 3D space.
	 * @see java.lang.Object#Object()
	 */
	public SimulationVirialWater(int nMolecules, double simTemperature) {
		this(new Space3D(), nMolecules, simTemperature);
	}
	
	/**
	 * Constructor for SimulationVirial.
	 */
	public SimulationVirialWater(Space space, int nMolecules, double simTemperature) {
		super(space);

		Default.TRUNCATE_POTENTIALS = false;
		
		final PhaseCluster phase = new PhaseCluster(this);
		phase.setBoundary(this.space.makeBoundary(Boundary.NONE));	
		species = new etomica.models.water.SpeciesWater(this);
		species.setNMolecules(nMolecules);
		elementCoordinator.go();
		
		Controller controller = new Controller(this);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(this);
		integrator = new IntegratorMC(this);
		integrator.setSleepPeriod(1);
		integrator.setDoSleep(true);
		integrator.setTemperature(simTemperature);
		integrator.setInterval(1);
//		MCMoveAtom mcMoveAtom = new MeterVirial.MyMCMoveAtom(integrator);
//		mcMoveAtom.setAdjustInterval(100000);
//		mcMoveAtom.setStepSize(2.0);
//		mcMoveAtom.setTunable(false);
		for(int n=1; n<nMolecules; n++) {
			MCMoveMoleculeMulti move = new MCMoveMoleculeMulti(integrator, n);
//			move.setAcceptanceTarget(0.10);
			move.setAdjustInterval(100000);
			move.setStepSize(2.0);
			move.setTunable(false);
		}
		MCMoveRotateMolecule3D moveRotate = new MCMoveRotateMolecule3D(integrator);
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(((AtomFactoryWater)species.moleculeFactory()).oFactory.type(),java.awt.Color.red);
		ColorSchemeByType.setColor(((AtomFactoryWater)species.moleculeFactory()).hFactory.type(),java.awt.Color.white);
		
		this.elementCoordinator.go();		
		
		Vector3D origin = new Vector3D(15.,15.,15.);
		SpeciesAgent speciesAgent = phase.getAgent(species);
		speciesAgent.coord.translateTo(new Vector3D(15.,15.,15.));
		speciesAgent.firstMolecule().coord.translateTo(new Vector3D(15.,15.,15.));
				
		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		AtomIteratorList list1 = new AtomIteratorList(childList);
		list1.reset();
		while(list1.hasNext()) list1.next().coord.translateTo(origin);
	}
	
	private MeterVirialAbs meterVirial;
	private double simTemperature;
	private double refTemperature;
	private PairSet pairs;
	private P0Cluster p2;
	private Species species;
	protected IntegratorMC integrator;
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterVirialAbs getMeterVirial() {
		return meterVirial;
	}

	/**
	 * Returns the refTemperature.
	 * @return double
	 */
	public double getRefTemperature() {
		return refTemperature;
	}

	/**
	 * Returns the simTemperature.
	 * @return double
	 */
	public double getSimTemperature() {
		return simTemperature;
	}

	/**
	 * Sets the meterVirial and configures displays for it.
	 * @param meterVirial The meterVirial to set
	 */
	public void setMeterVirialAbs(MeterVirialAbs meterVirial) {
		this.meterVirial = meterVirial;

		MeterDatumSourceWrapper bMeter = new MeterDatumSourceWrapper(meterVirial);
		bMeter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(this);
		bPlot.setDataSources(bMeter.getHistory());
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		bMeter.getHistory().setHistoryLength(1000);
		bPlot.setLabel("B running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(this);
		meterVirial.setHistorying(true);
		MeterDatumSourceWrapper[] clusterMeter = new MeterDatumSourceWrapper[4];
		for(int i=0; i<2; i++) {
			clusterMeter[i] = new MeterDatumSourceWrapper(meterVirial.allMeters()[i]);
			clusterMeter[i].setLabel("Meter"+i);
			clusterMeter[i].setHistorying(true);
			clusterMeter[i].getHistory().setHistoryLength(1000);
		}
		clusterPlot.setDataSources(new DataSource[] {
									clusterMeter[0].getHistory()//,
//									clusterMeter[1].getHistory()
									});
//		clusterPlot.setDataSources(meterVirial.getDataSources("History"));
//		clusterPlot.setWhichValue(MeterAbstract.AVERAGE);
		clusterPlot.setLabel("Cluster integrals");
		
		this.elementCoordinator.go();
		this.makeAndDisplayFrame();
	}

	/**
	 * Sets the refTemperature.
	 * @param refTemperature The refTemperature to set
	 */
	public void setRefTemperature(double refTemperature) {
		this.refTemperature = refTemperature;
	}

	/**
	 * Sets the simTemperature.
	 * @param simTemperature The simTemperature to set
	 */
	public void setSimTemperature(double simTemperature) {
		this.simTemperature = simTemperature;
	}
	
	public PairSet pairs() {
		return pairs;
	}
		
	public P0Cluster getSimPotential() {
		return p2;
	}
	
	public void setSimPotential(P0Cluster p2) {
		this.p2 = p2;
	}
	
	public Species species() {
		return species;
	}
		
	public static void main(String[] args) {
		double sigmaW = 3.167;
		double epsW = 78.21;
		double temperature = Kelvin.UNIT.toSim(600.0);  //temperature of calculated virial coefficient
		double simTemperature = 1.0*temperature; //temperature governing sampling of configurations
		double sigmaHSRef = 1.5*sigmaW*SimulationVirial.sigmaLJ1B(epsW/simTemperature);  //diameter of reference HS system
//		double sigmaHSMod = sigmaW*SimulationVirial.sigmaLJ1B(epsW/simTemperature); //range in which modified-f for sampling will apply abs() function
		System.out.println((epsW/simTemperature)+" "+sigmaHSRef);
//		System.out.println((epsW/temperature)+" "+sigmaHSMod);
		int nMolecules = 4;
		SimulationVirialWater sim = new SimulationVirialWater(nMolecules, simTemperature);
		
		//set up simulation potential
		P2WaterSPCE p2 = new P2WaterSPCE(new Simulation());//give dummy simulation
		Simulation.instance = sim;
//		MayerModified f = new MayerModified(p2LJ, sigmaHSMod);
		MayerFunction f = new MayerGeneral(p2);
		ClusterAbstract simCluster = null;
		double refIntegral = 1.0;
		switch (nMolecules) {
			case 2: simCluster = new B2(f);
					refIntegral = Standard.B2HS(sigmaHSRef); 
					break;
			case 3: simCluster = new C3(f); 
					refIntegral = Standard.C3HS(sigmaHSRef);
					break;
			case 4: simCluster = new ClusterSum(1.0,new Cluster[] {new D4(f, true), new D5(f, true), new D6(f)}); 
					break;
			default: throw new RuntimeException("Not able to handle nMolecules = "+nMolecules);
		}
		double refTemperature = temperature;
		Cluster refCluster = null;
		if(nMolecules == 2) refCluster = new B2(new MayerHardSphere(sigmaHSRef));
		else refCluster = new etomica.virial.cluster.Ring(nMolecules,1.0,new MayerHardSphere(sigmaHSRef));

		P0Cluster p0 = new P0Cluster(sim.hamiltonian.potential, simCluster);
		sim.setSimPotential(p0);
		
		MeterVirialAbs meterVirial = new MeterVirialAbs(sim,refTemperature, refCluster, refIntegral, p0);
		
		meterVirial.setBlockSize(10000);
		sim.setMeterVirialAbs(meterVirial);

		DisplayTable meterClusterTable = new DisplayTable();
		meterClusterTable.setDatumSources(meterVirial);
		meterClusterTable.setLabel("Cluster");
		meterClusterTable.setWhichValues(
			new DataSource.ValueType[] {MeterAbstract.CURRENT, MeterAbstract.AVERAGE, MeterAbstract.ERROR});
		meterClusterTable.setUpdateInterval(1000);
		meterClusterTable.setPrecision(8);
		DataSourceCountSteps meterCycles = new DataSourceCountSteps(sim);
		DisplayBox displayCycles = new DisplayBox(sim);
		displayCycles.setPrecision(10);
		displayCycles.setUpdateInterval(100000);
		meterCycles.setUpdateInterval(1);
		displayCycles.setDatumSource(meterCycles);
		
		ModifierStepSize modStep = sim.new ModifierStepSize();
		DeviceToggleButton stepSizeAdjustButton = new DeviceToggleButton(sim, modStep, "Adjusting step", "Not adjusting step");
		
		sim.elementCoordinator.go();
		meterClusterTable.setDatumSources(meterVirial.allMeters());
		meterClusterTable.addDatumSources(meterVirial);

	}//end of main
	
	private class ModifierStepSize implements ModifierBoolean {
		
		boolean value = false;
		public void setBoolean(boolean value) {
			this.value = value;
			MCMove[] moves = integrator.getMCMoves();
			for(int i=0; i<moves.length; i++) moves[i].setTunable(value);
		}
		public boolean getBoolean() {
			return integrator.getMCMoves()[0].getTunable();
		}
	}
	
}

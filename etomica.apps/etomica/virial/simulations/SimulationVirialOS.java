package etomica.virial.simulations;

import etomica.*;
import etomica.graphics.*;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.overlap.*;

/**
 * @author kofke
 *
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 */
public class SimulationVirialOS extends SimulationGraphic {

	/**
	 * Default constructor, using a 3D space.
	 * @see java.lang.Object#Object()
	 */
	public SimulationVirialOS(int nMolecules, double simTemperature) {
		this(new Space3D(), nMolecules, simTemperature);
	}
	
	/**
	 * Constructor for SimulationVirial.
	 */
	public SimulationVirialOS(Space space, int nMolecules, double simTemperature) {
		super(space);

		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
		
		phase = new Phase(this);
		phase.setBoundary(this.space.makeBoundary(Space3D.Boundary.NONE));	
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(nMolecules);
		elementCoordinator.go();
		pairs = new PairSet(((AtomTreeNodeGroup)phase.getAgent(species).node).childList);
		
		Controller controller = new Controller(this);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(this);
		integrator = new IntegratorMC(this);
		integrator.setSleepPeriod(1);
		integrator.setDoSleep(false);
		integrator.setTemperature(simTemperature);
		MCMoveAtom mcMoveAtom1 = new MeterVirial.MyMCMoveAtom(integrator);
		MCMoveAtomMulti mcMoveAtom2 = new MCMoveAtomMulti(integrator,2);
		if(nMolecules > 3) new MCMoveAtomMulti(integrator,3);
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		
		this.elementCoordinator.go();
		
//		Space3D.Vector origin = new Space3D.Vector(5.,5.,5.);
//		SpeciesAgent speciesAgent = phase.getAgent(species);
//		speciesAgent.coord.translateTo(new Space3D.Vector(5.,5.,5.));
//		speciesAgent.firstMolecule().coord.translateTo(new Space3D.Vector(5.,5.,5.));
//				
//		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
//		AtomIteratorList list1 = new AtomIteratorList(childList);
//		list1.reset();
//		while(list1.hasNext()) list1.next().coord.translateTo(origin);
	}
	
	public static double sigmaLJ1 = 1.28412293285;//  ( 4 + 4 sqrt(1-Ln(2)) ) / Ln(4))^(1/6), which is where f(r) = 1 for LJ

	/**
	 * Returns the separation r for the LJ potential at which beta*u(r) = -Ln(2)
	 * (so that f(r) = 1)
	 * @param beta
	 * @return double
	 */
	public static double sigmaLJ1B(double beta) {
		double log2 = Math.log(2.0);
		if(beta <= log2) return Math.pow(2.0, 1./6.);
		else return Math.pow( 2.0*(beta + Math.sqrt(beta*(beta-log2)) )/log2, 1./6.);
	}
	private MeterOverlap meter;
	private double simTemperature;
	private double refTemperature;
	private PairSet pairs;
	private P2Cluster p2;
	private SpeciesSpheresMono species;
	protected IntegratorMC integrator;
	private Phase phase; 
	
	public Phase phase() {return phase;}
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterOverlap getMeter() {
		return meter;
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
	public void setMeter(MeterOverlap meter) {
		this.meter = meter;

		MeterDatumSourceWrapper bMeter = new MeterDatumSourceWrapper(meter);
		bMeter.setHistorying(true);
		bMeter.getHistory().setNValues(1000);
		DisplayPlot bPlot = new DisplayPlot(this);
		bPlot.setDataSources(bMeter.getHistory());
		bPlot.setLabel("x=0 running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(this);
		clusterPlot.setDataSources(meter);
//		clusterPlot.setDataSources(meterVirial.getDataSources("History"));
		clusterPlot.setWhichValue(MeterAbstract.AVERAGE);
		clusterPlot.setLabel("Cluster integrals");
		
		this.elementCoordinator.go();
		clusterPlot.setDataSources(meter);
		bPlot.setDataSources(bMeter.getHistory());
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
		
	public P2Cluster getSimPotential() {
		return p2;
	}
	
	public void setSimPotential(P2Cluster p2) {
		this.p2 = p2;
	}
	
	public SpeciesSpheresMono species() {
		return species;
	}
		
	public static void main(String[] args) {
		double simTemperature = 1.3; //temperature governing sampling of configurations
		double sigmaHSRef = 1.0*sigmaLJ1B(1.0/simTemperature);  //diameter of reference HS system
		double b0 = 2*Math.PI/3. * Math.pow(sigmaHSRef,3);
		System.out.println("sigmaHSRef: "+sigmaHSRef);
		System.out.println("b0: "+b0);
		System.out.println("B3HS: "+(-5./8.*b0*b0));
//		sigmaHSRef = 1.0;
//		int nMolecules = 3;
		int nMolecules = 4;
		boolean targetPositive = true;   //calculating positive or negative part of target integral?
		boolean simulatingTarget = true; //simulating target system or HS reference?
		
		SimulationVirialOS sim = new SimulationVirialOS(nMolecules, simTemperature);
		
		ConfigurationCluster configuration = new ConfigurationCluster(sim);
		configuration.setPhase(sim.phase());
		sim.species().setDiameter(sigmaHSRef);
		sim.integrator.setDoSleep(false);
		
		
		MayerHardSphere fHS = new MayerHardSphere(sigmaHSRef);
		Cluster refCluster = new etomica.virial.cluster.Ring(nMolecules, 1.0, fHS);
		refCluster = new etomica.virial.cluster.ReeHoover(4, 1.0, new Cluster.BondGroup(fHS, Standard.D4));
		refCluster = new etomica.virial.Cluster(4, 1.0, new Cluster.BondGroup(fHS, Standard.D4));
		refCluster.setPairSet(sim.pairs());
		boolean refPositive = !refCluster.hasOddBondNumber();

		//set up simulation potential
		P2ClusterSigned p2 = new P2ClusterSigned(sim.hamiltonian.potential, sim.pairs());
		p2.setSignPositive(simulatingTarget ? targetPositive : refPositive);
		p2.setTemperature(simTemperature);		
		P2LennardJones p2LJ = new P2LennardJones(p2);
		MayerGeneral f = new MayerGeneral(p2LJ);
		Cluster targetCluster = new etomica.virial.cluster.Ring(nMolecules, 1.0, f);
		targetCluster = new etomica.virial.cluster.ReeHoover(4, 1.0, new Cluster.BondGroup(f, Standard.D4));
		targetCluster.setPairSet(sim.pairs());


		Cluster simCluster = simulatingTarget ? targetCluster : refCluster;
		Cluster nonSimCluster = simulatingTarget ? refCluster : targetCluster;
		p2.setCluster(simCluster);				
		sim.setSimPotential(p2);

		configuration.setCluster(simCluster);
		configuration.setSignPositive(simulatingTarget ? targetPositive : refPositive);
		sim.phase().setConfiguration(configuration);
		
		MeterOverlap meter = new MeterOverlap(sim, simCluster, nonSimCluster);
		meter.setTemperature(simTemperature);
		meter.setSignPositive(simulatingTarget ? refPositive : targetPositive);
		meter.setActive(true);
		sim.setMeter(meter);
		
		sim.makeAndDisplayFrame();
//		sim.elementCoordinator.go();
	}//end of main
}

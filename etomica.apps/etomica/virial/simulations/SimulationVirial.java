package etomica.virial.simulations;

import etomica.*;
import etomica.graphics.*;
import etomica.virial.*;
import etomica.virial.cluster.*;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class SimulationVirial extends SimulationGraphic {

	/**
	 * Default constructor, using a 3D space.
	 * @see java.lang.Object#Object()
	 */
	public SimulationVirial(int nMolecules, double simTemperature) {
		this(new Space3D(), nMolecules, simTemperature);
	}
	
	/**
	 * Constructor for SimulationVirial.
	 */
	public SimulationVirial(Space space, int nMolecules, double simTemperature) {
		super(space);

		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
		
		final PhaseCluster phase = new PhaseCluster(this);
		phase.setBoundary(this.space.makeBoundary(Space3D.Boundary.NONE));	
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(nMolecules);
		elementCoordinator.go();
		
		Controller controller = new Controller(this);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(this);
		integrator = new IntegratorMC(this);
		integrator.setSleepPeriod(1);
		integrator.setDoSleep(false);
		integrator.setTemperature(simTemperature);
		MCMoveAtom mcMoveAtom = new MeterVirial.MyMCMoveAtom(integrator);
		for(int n=2; n<nMolecules; n++) {
			MCMoveAtomMulti move = new MCMoveAtomMulti(integrator, n);
//			move.setAcceptanceTarget(0.10);
		}
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		
		this.elementCoordinator.go();
		
		Space3D.Vector origin = new Space3D.Vector(5.,5.,5.);
		SpeciesAgent speciesAgent = phase.getAgent(species);
		speciesAgent.coord.translateTo(new Space3D.Vector(5.,5.,5.));
		speciesAgent.firstMolecule().coord.translateTo(new Space3D.Vector(5.,5.,5.));
				
		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		AtomIteratorList list1 = new AtomIteratorList(childList);
		list1.reset();
		while(list1.hasNext()) list1.next().coord.translateTo(origin);
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
	private MeterVirial meterVirial;
	private double simTemperature;
	private double refTemperature;
	private PairSet pairs;
	private P0Cluster p2;
	private SpeciesSpheresMono species;
	protected IntegratorMC integrator;
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterVirial getMeterVirial() {
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
	public void setMeterVirial(MeterVirial meterVirial) {
		this.meterVirial = meterVirial;

		MeterDatumSourceWrapper bMeter = new MeterDatumSourceWrapper(meterVirial);
		bMeter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(this);
		bPlot.setDataSources(bMeter.getHistory());
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		bMeter.getHistory().setNValues(1000);
		bPlot.setLabel("B running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(this);
		meterVirial.setHistorying(true);
		MeterDatumSourceWrapper[] clusterMeter = new MeterDatumSourceWrapper[4];
		for(int i=0; i<2; i++) {
			clusterMeter[i] = new MeterDatumSourceWrapper(meterVirial.allMeters()[i]);
			clusterMeter[i].setLabel("Meter"+i);
			clusterMeter[i].setHistorying(true);
			clusterMeter[i].getHistory().setNValues(1000);
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
	
	public SpeciesSpheresMono species() {
		return species;
	}
		
	public static void main(String[] args) {
		double temperature = 1.3;  //temperature of calculated virial coefficient
		double simTemperature = 1.0*temperature; //temperature governing sampling of configurations
		double sigmaHSRef = 1.0;//*sigmaLJ1B(1.0/simTemperature);  //diameter of reference HS system
		double sigmaHSMod = sigmaLJ1B(1.0/simTemperature); //range in which modified-f for sampling will apply abs() function
		System.out.println((1./simTemperature)+" "+sigmaHSRef);
		System.out.println((1./temperature)+" "+sigmaHSMod);
		int nMolecules = 2;
		SimulationVirial sim = new SimulationVirial(nMolecules, simTemperature);
		
		sim.species().setDiameter(sigmaHSRef);
		
		//set up simulation potential
		P2LennardJones p2LJ = new P2LennardJones(new Simulation());//give dummy simulation
		Simulation.instance = sim;
		MayerModified f = new MayerModified(p2LJ, sigmaHSMod);
//		MayerFunction f = new etomica.virial.dos.MayerDOS2();
//		MayerFunction f = new MayerGeneral(p2LJ);
//		MayerE e = new MayerE(p2LJ);
//		MayerFunction f = new MayerHardSphere();
//		MayerE e = new MayerEHardSphere();
//		Cluster simCluster = new Ring(nMolecules, 1.0, f);
		Cluster simCluster = new Chain(nMolecules, 1.0, f);
//		Cluster simCluster = new C3e(f);
//		Cluster simCluster = new etomica.virial.dos.ClusterGrouped(7.0);
//		ClusterAbstract simCluster = new etomica.virial.cluster.ClusterBz(nMolecules, e);
//		ClusterAbstract simCluster = new etomica.virial.ClusterSum(1.0, Standard.B6Clusters(f));

		P0Cluster p0 = new P0Cluster(sim.hamiltonian.potential, simCluster);
		sim.setSimPotential(p0);
		
//		MeterVirial meterVirial = (nMolecules==3) ? (MeterVirial)new MeterVirialB3(sim, sigmaHSRef, p0, p2LJ)
//												  : (MeterVirial)new MeterVirialB4(sim, sigmaHSRef, p0, p2LJ);
		double refTemperature = temperature;
//		Cluster refCluster = new C3e(new MayerHardSphere(1.187193));
//		Cluster refCluster = new C3e(new MayerHardSphere(7.0));
//		ClusterAbstract refCluster = simCluster;
		Cluster refCluster = new etomica.virial.cluster.Full(nMolecules,1.0,new MayerHardSphere(1.0));
		System.out.println("Setting up target cluster");
//		ClusterAbstract targetCluster = new Cluster(nMolecules, 1.0, new Cluster.BondGroup[] {new Cluster.BondGroup(new MayerHardSphere(), Standard.chain(nMolecules))}, true);
		ClusterAbstract targetCluster = new Cluster(nMolecules, 1.0, new Cluster.BondGroup[] {new Cluster.BondGroup(new MayerGeneral(p2LJ), Standard.chain(nMolecules))}, true);
//		ClusterAbstract targetCluster = new Cluster(nMolecules, 1.0, new Cluster.BondGroup[] {new Cluster.BondGroup(f, Standard.chain(nMolecules))}, true);
//		ClusterAbstract targetCluster = new Cluster(nMolecules, 1.0, new Cluster.BondGroup[] {new Cluster.BondGroup(f, Standard.ring(nMolecules))}, true);
//		Cluster targetCluster = new C3e(f);
		double refIntegral = 1.0;
		MeterVirial meterVirial = new MeterVirial(sim, 
			refTemperature, refCluster, refIntegral, 
			new ClusterAbstract[] {targetCluster}, p0);
 		System.out.println("Done setting up meter");
//		System.out.println("TargetCluster\n"+targetCluster.toString());
//		Cluster dosCluster = new Ring(nMolecules, 1.0, new MayerE(p2LJ));
//		etomica.virial.dos.MeterDOS meterDOS = new etomica.virial.dos.MeterDOS(sim, dosCluster);
//		meterDOS.setP0(p0);
//		meterDOS.setLabel("DOS");
//		DisplayPlot histogramPlot = new DisplayPlot(sim);
//		histogramPlot.setDataSources(meterDOS);
//		histogramPlot.setLabel("DOS");
//		histogramPlot.setWhichValue(MeterAbstract.AVERAGE);
		
		meterVirial.setTemperature(temperature);
		meterVirial.setBlockSize(10000);
		sim.setMeterVirial(meterVirial);

		DisplayTable meterClusterTable = new DisplayTable();
		meterClusterTable.setDatumSources(meterVirial);
		meterClusterTable.setLabel("Cluster");
		meterClusterTable.setWhichValues(
			new DataSource.ValueType[] {MeterAbstract.CURRENT, MeterAbstract.AVERAGE, MeterAbstract.ERROR});
		meterClusterTable.setUpdateInterval(10000);
		meterClusterTable.setPrecision(8);
		//Evaluates virial using HS-excess cluster
		//Gamma = GammaHS(0) + (GammaHS(1)-GammaHS(0))*[<g-gHS0>_abs(g-gHS0)]/[<gHS1-gHS0>_abs(g-gHS0)]
		//Calculated here is Gamma - GammaHS(0)
//		double sigmaHSRef = sigmaLJ1B(1.0/simTemperature);  //diameter of reference HS system
//		int nMolecules = 3;
//		SimulationVirial sim = new SimulationVirial(nMolecules, simTemperature);
//		
//		sim.species().setDiameter(sigmaHSRef);
//		
//		//set up simulation potential
//		P2LennardJones p2LJ = new P2LennardJones(new Simulation());//give dummy simulation
//		Simulation.instance = sim;
//		MayerGeneral f = new MayerGeneral(p2LJ);
//		Cluster ring = new Ring(nMolecules, 1.0, f);
//		Cluster hs1 = new Cluster(new MayerHardSphere(sigmaHSRef), ring);
//		
//		Cluster simCluster = new ExcessHS(ring);
//		Cluster targetCluster = new ExcessHS(ring);
//		Cluster refCluster = new ExcessHS(hs1);
//		double refIntegral = Standard.C3HS(sigmaHSRef) - Standard.C3HS(1.0);
//		P0Cluster p0 = new P0Cluster(sim.hamiltonian.potential, simCluster);
//		sim.setSimPotential(p0);
//		
//		MeterVirial meterVirial = new MeterVirial(sim, 1.0, refCluster,
//					refIntegral, new Cluster[] {targetCluster}, p0);
//		meterVirial.setTemperature(temperature);
//		sim.setMeterVirial(meterVirial);
		sim.elementCoordinator.go();
		meterClusterTable.setDatumSources(meterVirial.allMeters());
		meterClusterTable.addDatumSources(meterVirial);
//		histogramPlot.setDataSources(meterDOS);


	}//end of main
}

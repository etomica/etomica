package etomica.virial;

import etomica.*;
import etomica.virial.cluster.*;
import etomica.units.Dimension;
import etomica.graphics.*;

/**
 * @author kofke
 *
 * Meter to compute the 4th virial coefficient.
 * allMeters[0] is the sum of the cluster diagrams, while allMeters[1], [2], [3]
 * correspond to the three cluster diagrams.
 * 
 */
public class MeterVirialB4 extends MeterVirial {
	
	/**
	 * Constructor for MeterVirialB4.
	 * @param sim
	 * @param pairSet
	 * @param refIntegral
	 * @param refCluster
	 * @param refMayer
	 * @param clusters
	 * @param mayer
	 * @param simulationPotential
	 */
	public MeterVirialB4(
		Simulation sim,
		PairSet pairSet,
		double refIntegral,
		Cluster refCluster,
		MayerFunction refMayer,
		Cluster[] clusters,
		MayerFunction mayer,
		P2Cluster simulationPotential) {
		super(
			sim,
			pairSet,
			refIntegral,
			refCluster,
			refMayer,
			clusters,
			mayer,
			simulationPotential);
	}

	/**
	 * Constructor for MeterVirialB4.
	 * @param sim
	 */
	
	public MeterVirialB4(Simulation sim) {
		super(sim, 4);
		clusters = new Cluster[3];
		clusters[0] = new D4();
		clusters[1] = new D5();
		clusters[2] = new D6();
		D5a = new Cluster(-3./4., new int[][] {{0,1},{0,3},{1,2},{1,3},{2,3}});
		setSigma(Default.ATOM_SIZE);
	}
	

	public static void main(String[] args) {
	
		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
//		Default.EXPLICIT_LOOP = true;
		
		SimulationGraphic sim = new SimulationGraphic(new Space3D());
		final Phase phase = new Phase(sim);
		phase.setBoundary(sim.space.makeBoundary(Space3D.Boundary.NONE));	
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		species.setNMolecules(4);
		Controller controller = new Controller(sim);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(sim);
		IntegratorMC integrator = new IntegratorMC(sim);
		integrator.setSleepPeriod(1);
		MyMCMoveAtom mcMoveAtom = new MyMCMoveAtom(integrator);
		
		P2MayerModified p2Mayer = new P2MayerModified(sim.hamiltonian.potential);
		P2LennardJones p2LJ = new P2LennardJones(p2Mayer);
		p2Mayer.setPotential(p2LJ);
		p2Mayer.setSpecies(species, species);
		
		MeterVirialB4 meterVirial = new MeterVirialB4(sim);
		
		double sigmaHS = 1.28412293285;  //  ( 4 + 4 sqrt(1-Ln(2)) ) / Ln(4))^(1/6), which is where f(r) = 1 for LJ
		double temperature = 1.3;
		
		p2Mayer.setSigma(sigmaHS);
		species.setDiameter(sigmaHS);
		
		p2Mayer.setTemperature(temperature);
		integrator.setTemperature(temperature);
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		

		AtomIteratorNeighbor nbrIterator1 = new AtomIteratorNeighbor(true);
		sim.elementCoordinator.go();
		
		Space3D.Vector origin = new Space3D.Vector(5.,5.,5.);
		SpeciesAgent speciesAgent = phase.getAgent(species);
		speciesAgent.coord.translateTo(new Space3D.Vector(5.,5.,5.));
		speciesAgent.firstMolecule().coord.translateTo(new Space3D.Vector(5.,5.,5.));
		
		meterVirial.setSigma(sigmaHS);
		meterVirial.setAtoms(((AtomTreeNodeGroup)speciesAgent.node).childList);
		meterVirial.setTemperature(temperature);
		meterVirial.setSimulationPotential(p2Mayer);
		
		NeighborManager.Criterion criterion = new NeighborManager.Criterion() {
			public boolean areNeighbors(Atom a1, Atom a2) {
				return Math.abs(a1.node.index()-a2.node.index()) == 1 ||
				   (a1==phase.firstAtom() && a2==phase.lastAtom()) ||
				   (a2==phase.firstAtom() && a1==phase.lastAtom());
			}};
//				return Math.abs(a1.node.index()-a2.node.index()) == 1;
//			}};
		
		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		nbrIterator1.setupNeighbors(childList, criterion);
		AtomIteratorList list1 = new AtomIteratorList(childList);
		AtomPairIterator api1= new ApiGeneral(sim.space,list1,nbrIterator1);
		p2Mayer.setIterator(api1);
		list1.reset();
		while(list1.hasNext()) list1.next().coord.translateTo(origin);
//		IteratorDirective.testSuitePair(api1,phase.firstAtom(),phase.firstAtom(),phase.lastAtom());

		MeterDatumSourceWrapper b4Meter = new MeterDatumSourceWrapper(meterVirial);
		b4Meter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(sim);
		bPlot.setDataSources(b4Meter.getHistory());
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		b4Meter.getHistory().setNValues(1000);
		bPlot.setLabel("B4 running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(sim);
		meterVirial.setHistorying(true);
		MeterDatumSourceWrapper[] clusterMeter = new MeterDatumSourceWrapper[4];
		for(int i=0; i<4; i++) {
			clusterMeter[i] = new MeterDatumSourceWrapper(meterVirial.allMeters()[i]);
			clusterMeter[i].setHistorying(true);
			clusterMeter[i].getHistory().setNValues(1000);
		}
		clusterMeter[0].setLabel("Reference");
		clusterPlot.setDataSources(new DataSource[] {
									clusterMeter[0].getHistory(),
									clusterMeter[1].getHistory(),
									clusterMeter[2].getHistory(),
									clusterMeter[3].getHistory()
									});
//		clusterPlot.setDataSources(meterVirial.getDataSources("History"));
//		clusterPlot.setWhichValue(MeterAbstract.AVERAGE);
		clusterPlot.setLabel("Cluster integrals");
		
		sim.elementCoordinator.go();

		sim.makeAndDisplayFrame();		
	
	}
}

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
public class MeterVirialB4 extends MeterGroup implements DatumSource {

	private Atom[] atoms;
	public static final int N = 4;
	private double sigma, sigma2;
	private double[][] f, fSim, fRef;
	private AtomPair[][] pairs;
	private Cluster[] clusters;
	private P2MayerModified simulationPotential;
	private double temperature, beta;
	private double B4Ref;
	private Cluster D5a;
	
	/**
	 * Constructor for MeterVirialB4.
	 * @param sim
	 */
	public MeterVirialB4(Simulation sim) {
		super(sim, 4);
		atoms = new Atom[N];
		f = new double[N-1][];
		fSim = new double[N-1][];
		fRef = new double[N-1][];
		pairs = new AtomPair[N-1][];
		for(int i=0; i<N-1; i++) {
			f[i] = new double[N-1-i];
			fSim[i] = new double[N-1-i];
			fRef[i] = new double[N-1-i];
			pairs[i] = new AtomPair[N-1-i];
			for(int j=0; j<N-1-i; j++) pairs[i][j] = new AtomPair(sim.space);
		}
//		clusters = new Cluster[3];
//		clusters[0] = new D4();
//		clusters[1] = new D5();
//		clusters[2] = new D6();
		clusters = new Cluster[3];
		clusters[0] = new D4();
		clusters[1] = new D5();
		clusters[2] = new D6();
		D5a = new Cluster(-3./4., new int[][] {{0,1},{0,3},{1,2},{1,3},{2,3}});
		setSigma(Default.ATOM_SIZE);
	}
	
	public void setAtoms(AtomList list) {
		AtomIteratorList iterator = new AtomIteratorList(list);
		iterator.reset();
		int i=0;
		while(iterator.hasNext() && i<N) atoms[i++] = iterator.next();
		if(atoms[N-1]==null) throw new RuntimeException("Not enough atoms given to MeterVirial");
		for(i=0; i<N-1; i++) {
			for(int j=0; j<N-1-i; j++) pairs[i][j].reset(atoms[i],atoms[i+j+1]);
		}
	}
	
	/**
	 * @see etomica.MeterGroup#updateValues()
	 */
	public void updateValues() {
		for(int i=0; i<N-1; i++) {
			for(int j=0; j<N-1-i; j++) {
				pairs[i][j].reset();
				fSim[i][j] = Math.exp(-beta*simulationPotential.energy(pairs[i][j]));
				double bu = simulationPotential.mostRecentBetaU();
				f[i][j] = Math.exp(-bu) - 1.0;
				fRef[i][j] = (pairs[i][j].r2()<sigma2) ? -1.0 : 0.0;
			}
		}
//		double DSim = fSim[0][0]*fSim[1][0]*fSim[2][0];  //chain
		double DSim = clusters[0].value(fSim);  //ring
		double D4 = clusters[0].value(f) - clusters[0].value(fRef);
		double D5 = 0.5*(clusters[1].value(f)+D5a.value(f)) - 0.5*(clusters[1].value(fRef)+D5a.value(fRef));
		double D6 = clusters[2].value(f) - clusters[2].value(fRef);
		double D4Ref = clusters[0].value(fRef);
		currentValues[3] = D4Ref/DSim;
		currentValues[0] = D4/DSim;
		currentValues[1] = D5/DSim;
		currentValues[2] = D6/DSim;
	}
	
	public double value(DataSource.ValueType type) {
		MeterScalar[] m = allMeters();
		double sum = 0.0;
		for(int i=0; i<clusters.length; i++) {
			sum += clusters[i].weight()*m[i].average();
		}
		double denom = m[3].average();
		return (denom != 0.0) ? B4Ref/m[3].average()*sum : Double.NaN;
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return Dimension.NULL;
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
	/**
	 * Returns the sigma.
	 * @return double
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Sets the sigma.
	 * @param sigma The sigma to set
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
		sigma2 = sigma*sigma;
		double b0 = 2.0*Math.PI/3.0*sigma2*sigma;
		B4Ref = 0.3238*8*b0*b0*b0;
	}
	
	/**
	 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
	 */
	private static class MyMCMoveAtom extends MCMoveAtom {
		MyMCMoveAtom(IntegratorMC integrator) {super(integrator);}
		
		public boolean doTrial() {
			if(phase.atomCount()==0) return false;
			atom = null;
			while(atom == null || atom.node.index()==0) atom = phase.speciesMaster.atomList.getRandom();
			uOld = potential.calculate(phase, iteratorDirective.set(atom), energy.reset()).sum();
			atom.coord.displaceWithin(stepSize);
			uNew = Double.NaN;
			return true;
		}
	}

	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
		beta = 1.0/temperature;
	}

	/**
	 * Returns the simulationPotential.
	 * @return P2MayerModified
	 */
	public P2MayerModified getSimulationPotential() {
		return simulationPotential;
	}

	/**
	 * Sets the simulationPotential.
	 * @param simulationPotential The simulationPotential to set
	 */
	public void setSimulationPotential(P2MayerModified simulationPotential) {
		this.simulationPotential = simulationPotential;
	}

}

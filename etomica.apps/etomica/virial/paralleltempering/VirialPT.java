package etomica.virial.paralleltempering;

import etomica.*;
import etomica.data.DataSourceAcceptanceRatio;
import etomica.data.DataSourceCountSteps;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPT;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.P2LennardJones;
import etomica.space3d.Boundary;
import etomica.space3d.Space3D;
import etomica.virial.*;
import etomica.virial.cluster.*;

public class VirialPT extends SimulationGraphic {
    
	PhaseCluster[] phase;
	SpeciesSpheresMono species;
	Controller controller;
	IntegratorMC[] integrator;
	P0Cluster potential;
    
    //computes temperature values in fixed ratios given max, min and number of values
	private static double[] tFixedRatio(double tMin, double tMax, int nT) {
		double[] temperatures = new double[nT];
		for(int i=0; i<nT; i++) {
			if(nT > 1) temperatures[i] = Math.exp(Math.log(tMin) + (double)i/(double)(nT-1)*Math.log(tMax/tMin));
			else temperatures[i] = tMin;
		}
		return temperatures;
	}

    
	public VirialPT(int nMolecules, double tMin, double tMax, int nPhase) {
		this(nMolecules, tFixedRatio(tMin, tMax, nPhase));
	}
	
	public VirialPT(int nMolecules, double[] temperatures) {
		super(new Space3D());
		Simulation.instance = this;
		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
		Default.HISTORY_PERIOD = 1000; //for phase tracking
        int nPhase = temperatures.length;
 
		P2LennardJones p2LJ = new P2LennardJones(new Simulation());//dummy simulation; parent of this potential will not be connected to the simulation
		Simulation.instance = this;
       
		double sigmaHSRef = 1.0;//etomica.virial.simulations.SimulationVirial.sigmaLJ1B(targetTemperature);
		double refTemperature = 1.0;
		double refIntegral = Standard.C3HS(sigmaHSRef);
//		Cluster[] clusters = MeterVirialB3.B3Clusters(p2LJ);
		Cluster[] clusters = MeterVirialB4.B4Clusters(p2LJ);
		Cluster refCluster = new Cluster(new MayerHardSphere(sigmaHSRef), clusters[0]);     
        
		elementCoordinator.addMediatorPair(new MediatorGraphic.DisplayMeter.NoAction(elementCoordinator));
        
		species = new SpeciesSpheresMono();
		species.setNMolecules(nMolecules);
		
		
		//take care to define other clusters (below) appropriately if using ReeHoover
		Cluster simCluster = clusters[0];
//		Cluster simCluster = new Ring(nMolecules, 1.0, f);
//		Cluster simCluster = new Cluster(5,1.0,
//			new Cluster.BondGroup(f, new int[][] {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{2,3},{2,4},{3,4}}));
//		Cluster simCluster = new ReeHoover(4, 1.0, new Cluster.BondGroup(f, Standard.D4));

//		P0Cluster potential = new P0Cluster(this.hamiltonian.potential);
//		potential.setCluster(simCluster);
		P0ClusterUmbrella potential = new P0ClusterUmbrella(this.hamiltonian.potential);
		Cluster[] umbClusters = new Cluster[clusters.length+1];
		umbClusters[0] = refCluster;
		for(int i=1; i<umbClusters.length; i++) umbClusters[i] = clusters[i-1];
		potential.setCluster(umbClusters);

		IntegratorPT.MCMoveSwapFactory swapFactory = new MCMoveSwapCluster.Factory();
		IntegratorPT integratorPT = new IntegratorPT(this, swapFactory);
		integratorPT.setDoSleep(true);
		integratorPT.setSleepPeriod(1);
		integratorPT.setInterval(species.getNMolecules());
		integratorPT.setSwapInterval(100);
        
		phase = new PhaseCluster[nPhase];

		integrator = new IntegratorMC[nPhase];
		DataSourceAcceptanceRatio[] meterAccept = new DataSourceAcceptanceRatio[nPhase-1];
		MeterVirial[] meterVirial = new MeterVirial[nPhase];
		double tLast = 0.0;
		for(int i=0; i<nPhase; i++) {
			double temperature = temperatures[i];
           
			phase[i] = new PhaseCluster();
			integrator[i] = new IntegratorMC(integratorPT);
			integrator[i].setTemperature(temperature);
			phase[i].setIntegrator(integrator[i]);
			phase[i].setBoundary(space.makeBoundary(Boundary.NONE));	
			elementCoordinator.go();
			
			ConfigurationCluster configuration = new ConfigurationCluster(this);
			configuration.setPhase(phase[i]);
			configuration.setCluster(simCluster);
			configuration.setSignPositive(true);
			phase[i].setConfiguration(configuration);						
            
			MCMoveAtom mcMoveAtom1 = new MeterVirial.MyMCMoveAtom(integrator[i]);
			for(int n=2; n<nMolecules; n++) {
				new MCMoveAtomMulti(integrator[i], n);
			}
			
			meterVirial[i] = new MeterVirial(this, 
				refTemperature, refCluster, refIntegral, 
				clusters, 
				potential);
            meterVirial[i].setPhase(phase[i]);
            meterVirial[i].setTemperature(temperature);
			meterVirial[i].setLabel("Temperature "+temperature);

			DisplayTable meterClusterTable = new DisplayTable();
			meterClusterTable.setDatumSources(meterVirial[i].allMeters());
			meterClusterTable.addDatumSources(meterVirial[i]);
			meterClusterTable.setLabel("Table "+temperature);
			meterClusterTable.setWhichValues(new DataSource.ValueType[] {MeterAbstract.CURRENT, MeterAbstract.AVERAGE, MeterAbstract.ERROR});
			meterClusterTable.setUpdateInterval(1000);
			meterClusterTable.setPrecision(5);
			integratorPT.addIntervalListener(meterClusterTable);
          
			integratorPT.addIntegrator(integrator[i]);
 
			DisplayPhase display = new DisplayPhase();
			display.setPhase(phase[i]);
			display.setLabel("Config "+Double.toString(integrator[i].temperature()));
            
			if(i>0) {
				meterAccept[i-1] = new DataSourceAcceptanceRatio(this, integratorPT.swapMoves()[i-1]);
				meterAccept[i-1].setLabel(Double.toString(tLast)+"<-->"+Double.toString(temperature));
			}
			tLast = temperature;
            
		}
        
		IntegratorPT.MeterPhaseTracker phaseTracker = integratorPT.new MeterPhaseTracker();
		phaseTracker.setPhases(phase);
		phaseTracker.setUpdateInterval(10);
		DisplayPlot phaseTrackPlot = new DisplayPlot();
		phaseTrackPlot.setDoLegend(false);
		phaseTrackPlot.setWhichValue(MeterAbstract.CURRENT);
		phaseTrackPlot.setDataSources(phaseTracker.dataSource());
		phaseTrackPlot.setLabel("Track");
		phaseTrackPlot.setUpdateInterval(10);
        
		DisplayTable acceptanceTable = new DisplayTable();
		acceptanceTable.setDatumSources(meterAccept);
		acceptanceTable.setLabel("Accept");
		acceptanceTable.setWhichValues(MeterAbstract.CURRENT);
        
//		DeviceButton logTableButton = new DeviceButton(perturbTable.makeLogTableAction());
//		logTableButton.setLabel("Perturb table out");
//        
//		DeviceButton logTableButton1 = new DeviceButton(m0Table.makeLogTableAction());
//		logTableButton1.setLabel("m0 table out");
                               
		controller = new Controller();
        
		DisplayBox cyclesBox = new DisplayBox(this, new DataSourceCountSteps());
		cyclesBox.setPrecision(6);
        
		ColorSchemeByType.setColor(species,java.awt.Color.green);
		DeviceTrioControllerButton buttons = new DeviceTrioControllerButton();
		
		int k = 0;
		for(k=0; k<meterVirial.length; k++) if(meterVirial[k].getTemperature()==1.3) break;
		MeterDatumSourceWrapper bMeter = new MeterDatumSourceWrapper(meterVirial[k]);
		bMeter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(this);
		bPlot.setDataSources(bMeter.getHistory());	
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		bMeter.getHistory().setHistoryLength(1000);
		bPlot.setLabel("B average ("+meterVirial[k].getTemperature()+")");
		
		DisplayTable virialTable = new DisplayTable(this);
		virialTable.setDatumSources(meterVirial);
		virialTable.setLabel("Virial");
		virialTable.setUpdateInterval(100);
		integratorPT.addIntervalListener(virialTable);
		
	}
        
	public static void main(String[] args) {
//		int nM = Integer.parseInt(args[0]);
//		double tMin = 1.0/Double.parseDouble(args[1]);
//		double tMax = 1.0/Double.parseDouble(args[2]);
//		int n = Integer.parseInt(args[3]);
//		VirialPT sim = new VirialPT(nM, tMin, tMax, n);
		int nM = 4;
		double[] temperature = new double[] {0.3, 0.7, 1.3, 2.5, 10.0};
		VirialPT sim = new VirialPT(nM, temperature);
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
	}
}
package etomica.association.GCPMWater;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.PNWaterGCPMReactionField;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.molecule.iterator.MoleculeIteratorAll;
import etomica.nbr.cell.molecule.BoxAgentSourceCellManagerMolecular;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.*;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.File;

/**
 * GCPM Water Monte Carlo NPT simulation in 3D.
 * average density = N*<1/V>
 * Translation/Rotation/Volume change
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 * @author Hye Min Kim
 * @author Andrew
 * @author navneeth
 */
public class SimGCPMWaterMCNVT extends Simulation {

    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public SpeciesGeneral species;
    public Box box;
    public PNWaterGCPMReactionField potential;



    public SimGCPMWaterMCNVT(int numMolceules, double densityMolLiter, double temperatureK, long numSteps) {
        super(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster();
        //setRandom(new RandomNumberGenerator(3));
        BoxAgentSourceCellManagerMolecular bASCellManagerMolecular = new BoxAgentSourceCellManagerMolecular(this, new MoleculePositionGeometricCenter(space), space);
        bASCellManagerMolecular.setRange(3.5);//association range=2.1-3.5
        System.out.println("numAtom=" +numMolceules);
        System.out.println("temperature = "+temperatureK+" K");
        System.out.println("numSteps = "+numSteps+" steps");
        System.out.println("density = "+densityMolLiter+" mol/L" +"\n");

    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        double density = rhoUnit.toSim(densityMolLiter);
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        double volume = 1/(density/numMolceules);
        double boxLength = Math.pow(volume, 1.0/3.0);      


	    mcMoveMolecule = new MCMoveMolecule(potentialMaster,random, space,10.0,15.0);
	    mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster,random, space);
	    box = new Box(space);
        addBox(box);
        integrator = new IntegratorMC(this.getRandom(), potentialMaster, box);
        integrator.setTemperature(temperature);

/*        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)mcMoveVolume.getTracker()).setNoisyAdjustment(true);*/

        integrator.getMoveManager().addMCMove(mcMoveMolecule);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        integrator.getMoveManager().setEquilibrating(true);
        this.getController().addActivity(new ActivityIntegrate(integrator), numSteps);
        //actionIntegrate.setSleepPeriod(1);
        species = SpeciesWater4P.create(new ConformationWaterGCPM(space));
        addSpecies(species);
        box.setNMolecules(species, numMolceules);
        BoxInflate inflater = new BoxInflate(box, space);//Performs actions that cause volume of system to expand or contract
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        
        potential = new PNWaterGCPMReactionField(space);//GCPMWater with Reaction Field Method
        potential.setBox(box);
        System.out.println("cut-off Distance "+boxLength*0.49+" A");
        potential.setTemperature(temperature);
        System.out.println("volume "+box.getBoundary().volume());
        System.out.println("Rho "+box.getMoleculeList().size() / box.getBoundary().volume()+"\n");
        
        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

        potentialMaster.addPotential(potential, new MoleculeIteratorAll(new ISpecies[]{species}, true), null);
        String configFile = "GCPM_NVT_"+numMolceules+"atoms"+temperatureK+"T"+numSteps+"steps";
        if(new File(configFile+".pos").exists()){
        	ConfigurationFile config = new ConfigurationFile(configFile);
            config.initializeCoordinates(box);
        }
        else{
	        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
	        config.initializeCoordinates(box);
        }
    }
 
    public static void main(String[] args) {  	

    	GCPMWaterMCParam params = new GCPMWaterMCParam();
    	ParseArgs.doParseArgs(params, args);
    	int numMolecules = params.numMolecules;
    	double density = params.density;
        double temperatureK = params.temperatureK;
        long numSteps = params.numSteps;
        double stepSizeTranslation = params.stepSizeTranslation;
        double stepSizeRotation = params.stepSizeRotation;

        long t1 = System.currentTimeMillis();

        final SimGCPMWaterMCNVT sim = new SimGCPMWaterMCNVT(numMolecules, density, temperatureK,numSteps);
        if (false) {
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE,"water", 1);
        	graphic.makeAndDisplayFrame();
        	return;
        }
        if (stepSizeTranslation==0.0){
        sim.integrator.getMoveManager().setEquilibrating(true);//adjusting step size automatically
        }
        else {
        	sim.integrator.getMoveManager().setEquilibrating(params.adjustStepSize);//if adjustStepSize is true the step size is adjusted automatically
        	sim.mcMoveMolecule.setStepSize(stepSizeTranslation);
        	sim.mcMoveRotateMolecule.setStepSize(stepSizeRotation);
        }

        final MeterPotentialEnergyFromIntegrator meterE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageCovariance energyAccumulator = new AccumulatorAverageCovariance(10);
        DataPumpListener energyManager = new DataPumpListener(meterE, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyManager);

        if (false) {
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE,"water", 1);
        	ISpecies species = sim.getSpecies(0);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
        	AccumulatorHistory densityHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        	DisplayPlot rhoPlot = new DisplayPlot();
        	densityHistory.setDataSink(rhoPlot.getDataSet().makeDataSink());
        	rhoPlot.setLabel("density");
        	graphic.add(rhoPlot);
        	DataSourceCountSteps stepCounter = new DataSourceCountSteps(sim.integrator);
        	AccumulatorHistory energyHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        	energyAccumulator.addDataSink(energyHistory, new StatType[]{energyAccumulator.MOST_RECENT});
        	DisplayPlot energyPlot = new DisplayPlot();
        	energyHistory.setTimeDataSource(stepCounter);
        	energyHistory.setDataSink(energyPlot.getDataSet().makeDataSink());
        	energyPlot.setLabel("energy");
        	graphic.add(energyPlot);
        	graphic.makeAndDisplayFrame();
        	return;
        }
        sim.integrator.reset();
        double initialEnthalpy = meterE.getData().getValue(0);
        System.out.println("initial Energy "+initialEnthalpy+"\n");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));
        System.out.println("step size of mcMoveMolecule "+sim.mcMoveMolecule.getStepSize());
        System.out.println("step size of mcMoveRotateMolecule "+sim.mcMoveRotateMolecule.getStepSize());
        
        System.out.println("Acceptance of mcMoveMolecule "+sim.mcMoveMolecule.getTracker().acceptanceProbability());
        System.out.println("Acceptance of mcMoveRotateMolecule "+sim.mcMoveRotateMolecule.getTracker().acceptanceProbability()+"\n");
        WriteConfiguration writeConfig = new WriteConfiguration(sim.space);
        writeConfig.setDoApplyPBC(false);//some atoms outside the box are allowed. Don't move atoms outside the box to the other side of the box.
        writeConfig.setBox(sim.box);
        writeConfig.setFileName("GCPM_NVT_"+numMolecules+"atoms"+temperatureK+"T"+numSteps+"steps.pos");
        writeConfig.actionPerformed();

        double finalEnergy = meterE.getData().getValue(0);
        System.out.println("final Energy "+finalEnergy+"\n");

        double avgEnergy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index).getValue(0);
        System.out.println("average Energy (simE) = "+avgEnergy);
        double varEnergy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.COVARIANCE.index).getValue(0);
        System.out.println("varEnergy = "+varEnergy+"\n");

        double temperature = Kelvin.UNIT.toSim(temperatureK);

        double avgCv = varEnergy/(temperature*temperature);
        System.out.println("average Cv(sim) = "+avgCv);
        System.out.println("average Cv(sim/n) = "+avgCv/numMolecules);
        Unit JoulesPerMolesK = new CompoundUnit(new Unit[]{Joule.UNIT,Mole.UNIT,Kelvin.UNIT},new double[]{1.0,-1.0,-1.0});
        double unitCvJpmK = JoulesPerMolesK.fromSim(1);
        System.out.println("unit conv Cv = "+unitCvJpmK+"\n");

        long t2 = System.currentTimeMillis();
        System.out.println("time = "+(t2-t1)/1000.0);

    }

    public static class GCPMWaterMCParam extends ParameterBase {
		public int numMolecules = 256;
		public double density = 1.00904;             //mol/L
		public double temperatureK = 600.0;          //Kelvin
		public long numSteps = 10000;
		public double stepSizeTranslation = 0.0;
		public double stepSizeRotation = 0.0;
		public boolean adjustStepSize = false;
	}

}

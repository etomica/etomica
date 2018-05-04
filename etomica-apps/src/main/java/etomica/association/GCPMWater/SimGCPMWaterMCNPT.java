package etomica.association.GCPMWater;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.PNWaterGCPMReactionField;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.nbr.cell.molecule.BoxAgentSourceCellManagerMolecular;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
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
public class SimGCPMWaterMCNPT extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public MCMoveVolume mcMoveVolume;
    public SpeciesWater4P species;
    public Box box;
    public PNWaterGCPMReactionField potential;
    double epsilon = 1.0;
    public ActivityIntegrate actionIntegrator;
        
    
    public SimGCPMWaterMCNPT(int numAtoms, double pressureBar, double densityMolLiter, double temperatureK, long numSteps) {
        super(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster();
        //setRandom(new RandomNumberGenerator(3));
        BoxAgentSourceCellManagerMolecular bASCellManagerMolecular = new BoxAgentSourceCellManagerMolecular(this, new MoleculePositionGeometricCenter(space), space);
        bASCellManagerMolecular.setRange(3.5);//association range=2.1-3.5
        BoxAgentManager cellAgentManager = new BoxAgentManager(bASCellManagerMolecular,NeighborCellManagerMolecular.class);
        System.out.println("pressure = "+pressureBar+"bar");
        System.out.println("initial density = "+densityMolLiter+"mol/L");
        System.out.println("temperature = "+temperatureK+"K");
        System.out.println("numSteps = "+numSteps+"steps");
    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        double density = rhoUnit.toSim(densityMolLiter);
        double pressure = Bar.UNIT.toSim(pressureBar);
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        double volume = 1/(density/numAtoms);
        double boxLength = Math.pow(volume, 1.0/3.0);      

	    integrator = new IntegratorMC(this, potentialMaster);
	    integrator.setTemperature(temperature);
	    mcMoveMolecule = new MCMoveMolecule(potentialMaster,random, space,10.0,15.0);
	    mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster,random, space);
	    mcMoveVolume = new MCMoveVolume(this,potentialMaster,space);//volume change
	    mcMoveVolume.setPressure(pressure);
	    box = new Box(space);
        addBox(box);
        Unit calPerMole = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});

        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)mcMoveVolume.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(mcMoveMolecule);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        integrator.getMoveManager().addMCMove(mcMoveVolume);
        integrator.getMoveManager().setEquilibrating(true);
        actionIntegrator = new ActivityIntegrate(integrator);
        //actionIntegrate.setSleepPeriod(1);
        actionIntegrator.setMaxSteps(numSteps);
        getController().addAction(actionIntegrator);
        species = new SpeciesWater4P(space);
        addSpecies(species);
        species.setConformation(new ConformationWaterGCPM(space));
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);//Performs actions that cause volume of system to expand or contract
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        
        potential = new PNWaterGCPMReactionField(space);//GCPMWater with Reaction Field Method
        potential.setBox(box);
        System.out.println("cut-off Distance "+boxLength*0.49+" A");
        potential.setTemperature(temperature);

        System.out.println("number of molecules "+box.getMoleculeList().getMoleculeCount());
        System.out.println("volume "+box.getBoundary().volume());
        System.out.println("Rho "+box.getMoleculeList().getMoleculeCount() / box.getBoundary().volume());
        
        BoxImposePbc pbc = new BoxImposePbc(box, space);
        pbc.setApplyToMolecules(true);
        integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

        potentialMaster.addPotential(potential, new ISpecies[] {species});
        String configFile = "GCPM_NPT_"+numAtoms+"atoms"+temperatureK+"T"+pressureBar+"Bar"+numSteps+"steps";
        if(new File(configFile+".pos").exists()){
        	ConfigurationFile config = new ConfigurationFile(configFile);
            config.initializeCoordinates(box);
        }
        else{
	        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
	        config.initializeCoordinates(box);
        }
        integrator.setBox(box);

    }
 
    public static void main(String[] args) {  	

    	GCPMWaterMCParam params = new GCPMWaterMCParam();
    	ParseArgs.doParseArgs(params, args);
    	int numAtoms = params.numAtoms;
    	double pressure = params.pressure;
    	double density = params.density;
        double temperature = params.temperature;
        long numSteps = params.numSteps;
        double stepSizeTranslation = params.stepSizeTranslation;
        double stepSizeRotation = params.stepSizeRotation;
        double stepSizeVolume = params.stepSizeVolume;

        final SimGCPMWaterMCNPT sim = new SimGCPMWaterMCNPT(numAtoms, pressure, density, temperature,numSteps);
        if (false) {
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE,"water", 1, sim.space,sim.getController());
        	graphic.makeAndDisplayFrame();
        	sim.actionIntegrator.setMaxSteps(Long.MAX_VALUE);
        	return;
        }
        if (stepSizeTranslation==0.0){
        sim.integrator.getMoveManager().setEquilibrating(true);//adjusting step size automatically
        }
        else {
        	sim.integrator.getMoveManager().setEquilibrating(params.adjustStepSize);//if adjustStepSize is true the step size is adjusted automatically
        	sim.mcMoveMolecule.setStepSize(stepSizeTranslation);
        	sim.mcMoveRotateMolecule.setStepSize(stepSizeRotation);
        	sim.mcMoveVolume.setStepSize(stepSizeVolume);
        }
        
        sim.actionIntegrator.setMaxSteps(numSteps);
        MeterDensity rhoMeter = new MeterDensity(sim.space);//Meter for measurement of the total molecule number density((number of molecules)/(volume of box)) in a box 
        rhoMeter.setBox(sim.box);
        AccumulatorAverage rhoAccumulator = new AccumulatorAverageFixed(1000);//Accumulator that keeps statistics for averaging and error analysis
        DataPump rhoPump = new DataPump(rhoMeter,rhoAccumulator);
        IntegratorListenerAction listener = new IntegratorListenerAction(rhoPump);
        listener.setInterval(1);
        sim.integrator.getEventManager().addListener(listener);
        final MeterEnthalpyVolume meterEV = new  MeterEnthalpyVolume(sim.integrator, pressure);
        AccumulatorAverageCovariance energyAccumulator = new AccumulatorAverageCovariance(10);
        DataPumpListener energyManager = new DataPumpListener(meterEV, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyManager);

        if (false) {
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE,"water", 1, sim.space,sim.getController());
        	SpeciesWater4P species = (SpeciesWater4P)sim.getSpecies(0);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
        	AccumulatorHistory densityHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        	rhoAccumulator.addDataSink(densityHistory, new StatType[]{rhoAccumulator.MOST_RECENT});
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
        	sim.actionIntegrator.setMaxSteps(Long.MAX_VALUE);
        	return;
        }
        double initialEnergy = meterEV.getData().getValue(0);
        System.out.println("initialEnergy "+initialEnergy);
        sim.getController().actionPerformed();
        System.out.println("step size of mcMoveMolecule "+sim.mcMoveMolecule.getStepSize());
        System.out.println("step size of mcMoveRotateMolecule "+sim.mcMoveRotateMolecule.getStepSize());
        System.out.println("step size of mcMoveVolume "+sim.mcMoveVolume.getStepSize());
        
        System.out.println("Acceptance of mcMoveMolecule "+sim.mcMoveMolecule.getTracker().acceptanceProbability());
        System.out.println("Acceptance of mcMoveRotateMolecule "+sim.mcMoveRotateMolecule.getTracker().acceptanceProbability());
        System.out.println("Acceptance of mcMoveVolume "+sim.mcMoveVolume.getTracker().acceptanceProbability());
        
        WriteConfiguration writeConfig = new WriteConfiguration(sim.space);
        writeConfig.setDoApplyPBC(false);//some atoms outside the box are allowed. Don't move atoms outside the box to the other side of the box.
        writeConfig.setBox(sim.box);
        writeConfig.setFileName("GCPM_NPT_"+numAtoms+"atoms"+temperature+"T"+pressure+"Bar"+numSteps+"steps.pos");
        writeConfig.actionPerformed();
        
    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        double finalDensity = rhoMeter.getDataAsScalar();
        finalDensity = rhoUnit.fromSim(finalDensity);
        System.out.println("next initial Density "+finalDensity);
        System.out.println("numAtom=" +numAtoms);

        double avgDensity = ((DataDouble)((DataGroup)rhoAccumulator.getData()).getData(rhoAccumulator.AVERAGE.index)).x;//average density
        double errDensity = ((DataDouble)((DataGroup)rhoAccumulator.getData()).getData(rhoAccumulator.ERROR.index)).x;
        double correlationBlock = ((DataDouble)((DataGroup)rhoAccumulator.getData()).getData(rhoAccumulator.BLOCK_CORRELATION.index)).x;

        System.out.println("err "+errDensity);
        System.out.println("correlationBlock "+correlationBlock);
        System.out.println("average density(sim)= " +avgDensity);
        avgDensity = rhoUnit.fromSim(avgDensity);
        System.out.println("average density(mol/liter)= " +avgDensity);

        double Z = pressure/(avgDensity*sim.integrator.getTemperature());
        double avgEnthalpy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index).getValue(0);
        System.out.println("average enthalpy (sim) = "+avgEnthalpy);
        double Enthalpy = meterEV.getData().getValue(0);
        System.out.println("final enthalpy = "+Enthalpy);
        avgEnthalpy /= numAtoms;
        System.out.println("Z="+Z);
    	Unit JoulesPerMoles = new CompoundUnit(new Unit[]{Joule.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
    	System.out.println("Average Enthalpy (Joule/mole) = "+JoulesPerMoles.fromSim(avgEnthalpy));
        double temp = sim.integrator.getTemperature();
        double varEnthalpy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.COVARIANCE.index).getValue(0);
        double varV = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.COVARIANCE.index).getValue(3);
        double covVEnthalpy = ((DataGroup)energyAccumulator.getData()).getData((energyAccumulator).COVARIANCE.index).getValue(1);

        System.out.println("varEnthalpy = "+varEnthalpy );
        System.out.println("varV = "+varV );
        System.out.println("covVEnthalpy = "+covVEnthalpy );

    }

    public static class GCPMWaterMCParam extends ParameterBase {
		public int numAtoms = 20;
		public double pressure = 170;//bar
		public double density = 3.3531383834004864;//1g/cm3=1000/18.02mol/L
		public double temperature = 630.0;//Kelvin
		public long numSteps = 5000;
		public double stepSizeTranslation = 0.0;
		public double stepSizeRotation = 0.0;
		public double stepSizeVolume = 0.0;
		public boolean adjustStepSize = false;
	}

}

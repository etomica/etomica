package etomica.association.GCPMWater;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.WriteConfiguration;

import etomica.action.activity.ActivityIntegrate2;
import etomica.box.Box;
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
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveVolume;
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

        
    
    public SimGCPMWaterMCNPT(int numMolceules, double pressureBar, double densityMolLiter, double temperatureK, long numSteps) {
        super(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster();
        //setRandom(new RandomNumberGenerator(3));
        BoxAgentSourceCellManagerMolecular bASCellManagerMolecular = new BoxAgentSourceCellManagerMolecular(this, new MoleculePositionGeometricCenter(space), space);
        bASCellManagerMolecular.setRange(3.5);//association range=2.1-3.5
        System.out.println("numAtom=" +numMolceules);
        System.out.println("temperature = "+temperatureK+" K");
        System.out.println("pressure = "+pressureBar+" bar");
        System.out.println("numSteps = "+numSteps+" steps");
        System.out.println("initial density = "+densityMolLiter+" mol/L" +"\n");

    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        double density = rhoUnit.toSim(densityMolLiter);
        double pressure = Bar.UNIT.toSim(pressureBar);
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        double volume = 1/(density/numMolceules);
        double boxLength = Math.pow(volume, 1.0/3.0);      


	    mcMoveMolecule = new MCMoveMolecule(potentialMaster,random, space,10.0,15.0);
	    mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster,random, space);
	    mcMoveVolume = new MCMoveVolume(this,potentialMaster,space);//volume change
	    mcMoveVolume.setPressure(pressure);
	    box = new Box(space);
        addBox(box);
        integrator = new IntegratorMC(this, potentialMaster,box);
        integrator.setTemperature(temperature);
        Unit calPerMole = new CompoundUnit(new Unit[]{Calorie.UNIT,Mole.UNIT},new double[]{1.0,-1.0});

/*        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)mcMoveVolume.getTracker()).setNoisyAdjustment(true);*/

        integrator.getMoveManager().addMCMove(mcMoveMolecule);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        integrator.getMoveManager().addMCMove(mcMoveVolume);
        integrator.getMoveManager().setEquilibrating(true);
        this.getController().addActivity(new ActivityIntegrate2(integrator), numSteps);
        //actionIntegrate.setSleepPeriod(1);
        species = new SpeciesWater4P(space);
        addSpecies(species);
        species.setConformation(new ConformationWaterGCPM(space));
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
        String configFile = "GCPM_NPT_"+numMolceules+"atoms"+temperatureK+"T"+pressureBar+"Bar"+numSteps+"steps";
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
    	double pressure = params.pressure;
    	double density = params.density;
        double temperature = params.temperature;
        long numSteps = params.numSteps;
        double stepSizeTranslation = params.stepSizeTranslation;
        double stepSizeRotation = params.stepSizeRotation;
        double stepSizeVolume = params.stepSizeVolume;

        long t1 = System.currentTimeMillis();

        final SimGCPMWaterMCNPT sim = new SimGCPMWaterMCNPT(numMolecules, pressure, density, temperature,numSteps);
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
        	sim.mcMoveVolume.setStepSize(stepSizeVolume);
        }

        MeterDensity rhoMeter = new MeterDensity(sim.box);
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
        	SimulationGraphic graphic = new SimulationGraphic(sim,SimulationGraphic.TABBED_PANE,"water", 1);
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
        	return;
        }
        sim.integrator.reset();
        double initialEnergy = sim.integrator.getPotentialEnergy();
        System.out.println("initial Energy "+initialEnergy);
        double initialEnthalpy = meterEV.getData().getValue(0);
        System.out.println("initial Enthaply "+initialEnthalpy+"\n");
sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);
        System.out.println("step size of mcMoveMolecule "+sim.mcMoveMolecule.getStepSize());
        System.out.println("step size of mcMoveRotateMolecule "+sim.mcMoveRotateMolecule.getStepSize());
        System.out.println("step size of mcMoveVolume "+sim.mcMoveVolume.getStepSize());
        
        System.out.println("Acceptance of mcMoveMolecule "+sim.mcMoveMolecule.getTracker().acceptanceProbability());
        System.out.println("Acceptance of mcMoveRotateMolecule "+sim.mcMoveRotateMolecule.getTracker().acceptanceProbability());
        System.out.println("Acceptance of mcMoveVolume "+sim.mcMoveVolume.getTracker().acceptanceProbability()+"\n");
        WriteConfiguration writeConfig = new WriteConfiguration(sim.space);
        writeConfig.setDoApplyPBC(false);//some atoms outside the box are allowed. Don't move atoms outside the box to the other side of the box.
        writeConfig.setBox(sim.box);
        writeConfig.setFileName("GCPM_NPT_"+numMolecules+"atoms"+temperature+"T"+pressure+"Bar"+numSteps+"steps.pos");
        writeConfig.actionPerformed();
        
    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        double finalDensity = rhoMeter.getDataAsScalar();
        finalDensity = rhoUnit.fromSim(finalDensity);
        System.out.println("next initial Density "+finalDensity+"\n");

        double avgDensity = ((DataDouble)((DataGroup)rhoAccumulator.getData()).getData(rhoAccumulator.AVERAGE.index)).x;//average density
        double errDensity = ((DataDouble)((DataGroup)rhoAccumulator.getData()).getData(rhoAccumulator.ERROR.index)).x;
        double correlationBlock = ((DataDouble)((DataGroup)rhoAccumulator.getData()).getData(rhoAccumulator.BLOCK_CORRELATION.index)).x;

        System.out.println("average density(sim)= " +avgDensity);
        System.out.println("err "+errDensity);
        System.out.println("correlationBlock "+correlationBlock);
        double avgDensitymolpl = rhoUnit.fromSim(avgDensity);
        System.out.println("average density(mol/liter)= " +avgDensitymolpl+"\n");

        double finalEnergy = sim.integrator.getPotentialEnergy();
        System.out.println("final Energy "+finalEnergy);
        double finalEnthalpy = meterEV.getData().getValue(0);
        System.out.println("final Enthalpy = "+finalEnthalpy+"\n");

        double Z = pressure/(avgDensity*sim.integrator.getTemperature());
        System.out.println("Z="+Z+"\n");

        double avgEnthalpy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index).getValue(0);
        System.out.println("average enthalpy (sim) = "+avgEnthalpy);
        avgEnthalpy /= numMolecules;
    	Unit JoulesPerMoles = new CompoundUnit(new Unit[]{Joule.UNIT,Mole.UNIT},new double[]{1.0,-1.0});
    	double avgEnthalpyjpm = JoulesPerMoles.fromSim(avgEnthalpy);
    	System.out.println("Average Enthalpy (Joule/mole) = "+avgEnthalpyjpm+"\n");

        double Volume = meterEV.getData().getValue(1);
        System.out.println("final volume = "+Volume);
        double avgVolume = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index).getValue(1);
        System.out.println("average volume (sim) = "+avgVolume);
        double avgVolumem3 = avgVolume*1e-27;
        System.out.println("Average Volume (Liter) = "+avgVolumem3+"\n");

        double temp = sim.integrator.getTemperature();
        System.out.println("Tsim = "+temp);
        double varEnthalpy = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.COVARIANCE.index).getValue(0);
        double varV = ((DataGroup)energyAccumulator.getData()).getData(energyAccumulator.COVARIANCE.index).getValue(3);
        double covVEnthalpy = ((DataGroup)energyAccumulator.getData()).getData((energyAccumulator).COVARIANCE.index).getValue(1);

        System.out.println("varEnthalpy = "+varEnthalpy );
        System.out.println("varV = "+varV );
        System.out.println("covVEnthalpy = "+covVEnthalpy+"\n");

        double avgKt = varV/(avgVolume*1*temp);
        double avgCp = varEnthalpy/(1*(temp*temp))/numMolecules;
        double avgalphap = covVEnthalpy/(1*(temp*temp)*avgVolume);
        double avgCv = avgCp - ((temp*(avgalphap*avgalphap))/(avgDensity*avgKt));

        System.out.println("average Kt(sim) = "+avgKt);
        System.out.println("average Cp(sim) = "+avgCp);
        System.out.println("average alphap(sim) = "+avgalphap+"\n");

        System.out.println("average Cv(sim) = "+avgCv);
        Unit JoulesPerMolesK = new CompoundUnit(new Unit[]{Joule.UNIT,Mole.UNIT,Kelvin.UNIT},new double[]{1.0,-1.0,-1.0});
        double avgCvJpmK = JoulesPerMolesK.fromSim(avgCv);
        System.out.println("average Cv(Joules/mol/K) = "+avgCvJpmK+"\n");
        System.out.println("conv Cv from sim to JpmK = "+JoulesPerMolesK.fromSim(1)+"\n");

        long t2 = System.currentTimeMillis();
        System.out.println("time = "+(t2-t1)/1000.0+"\n");

    }

    public static class GCPMWaterMCParam extends ParameterBase {
		public int numMolecules = 256;
		public double pressure = 123;//bar
		public double density = 3.3531383834004864;//1g/cm3=1000/18.02mol/L
		public double temperature = 800.0;//Kelvin
		public long numSteps = 10;
		public double stepSizeTranslation = 0.0;
		public double stepSizeRotation = 0.0;
		public double stepSizeVolume = 0.0;
		public boolean adjustStepSize = false;
	}

}

package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;

import etomica.action.BoxImposePbc;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.RandomNumberGenerator;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom on a surface, modeled with LJ.
 * 
 * @author msellers
 *
 */

public class SimDimerLJadatom extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerLJadatom";
    public final PotentialMaster potentialMaster;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorDimerMin;
    public IBox box;
    public IVector [] saddle, normal;
    public SpeciesSpheresMono fixed, movable;
    public P2LennardJones potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    public CalcGradientDifferentiable calcGradientDifferentiable;
    public CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    public IAtomSet movableSet;
    public Boolean saddleFine, calcModes, minSearch, normalDir;
    

    public SimDimerLJadatom(String fileName, Boolean useConfig, Boolean ortho, Boolean saddleFine, Boolean calcModes, Boolean minSearch, Boolean normalDir) {
    	super(Space3D.getInstance(), true);
    	potentialMaster = new PotentialMaster(space);
    	
    //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(random, 0, 5, space), space);
        addBox(box);
        
        BoxImposePbc imposePbc = new BoxImposePbc(box, space);
        
    // INTEGRATOR - MD
    	integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, space);
    	integratorMD.setTimeStep(0.01);
    	integratorMD.setTemperature(0.1);
    	integratorMD.setThermostatInterval(100);
    	integratorMD.setIsothermal(true);
    	integratorMD.setBox(box);
    	activityIntegrateMD = new ActivityIntegrate(integratorMD);
    	integratorMD.addIntervalAction(imposePbc);
    	
    //SPECIES
    	double sigma = 1.0;
    	Tin tinFixed = new Tin("SnFixed", Double.POSITIVE_INFINITY);
    	fixed = new SpeciesSpheresMono(this, space, tinFixed);
        movable = new SpeciesSpheresMono(this, space);      
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)fixed.getLeafType()).setDiameter(sigma);
        ((AtomTypeSphere)movable.getLeafType()).setDiameter(sigma);
    	
        // Must be in same order as the respective species is added to SpeciesManager
        box.setNMolecules(fixed, 256);    	
    	
    	box.setDensity(1);
    	
    	potential = new P2LennardJones(space, sigma, 1.0);
		potentialMaster.addPotential(potential, new AtomType[]{fixed.getLeafType(), fixed.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{movable.getLeafType(), fixed.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{movable.getLeafType(), movable.getLeafType()});
        
    //CRYSTAL
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        config.initializeCoordinates(box); 
       
    //ADATOM CREATION AND PLACEMENT
        IMolecule iMolecule = movable.makeMolecule();
        box.addMolecule(iMolecule);
        IVector adAtomPos = ((IAtomPositioned)iMolecule.getChildList().getAtom(0)).getPosition();
        adAtomPos.setX(0, box.getBoundary().getDimensions().x(0)/2+0.5);
        adAtomPos.setX(1, box.getBoundary().getDimensions().x(0)/16 + 0.1);
        adAtomPos.setX(2, box.getBoundary().getDimensions().x(0)/16);
        
  //INTEGRATOR - Dimer
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new ISpecies[]{movable}, space);
        integratorDimer.setBox(box);
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
        integratorDimer.setOrtho(ortho, false, false);
        integratorDimer.setFileName(fileName);
        integratorDimer.setActivityIntegrate(activityIntegrateDimer);

    //ADD CONTROLLER ACTIONS
    	getController().addAction(activityIntegrateMD);
    	getController().addAction(activityIntegrateDimer);
    	
    //FINE-DIMER SETTINGS
        if(saddleFine){
        	getController().removeAction(activityIntegrateMD);
        	
        	ConfigurationFile configFile = new ConfigurationFile(fileName+"_saddle");
        	configFile.initializeCoordinates(box);
        	
        	integratorDimer.file = fileName+"_fine";
            integratorDimer.deltaR = 0.0005;
            integratorDimer.dXl = 10E-5;       
            integratorDimer.deltaXmax = 0.005;
            integratorDimer.dFsq = 0.0001*0.0001;
            integratorDimer.dFrot = 0.01;
        }
                        
    //INTEGRATOR - Minimum Energy Path
        if(minSearch){
        	ConfigurationFile configFile = new ConfigurationFile(fileName+"_fine_saddle");
        	configFile.initializeCoordinates(box);
            integratorDimerMin = new IntegratorDimerMin(this, potentialMaster, new ISpecies[]{movable}, fileName, normalDir, space);
            integratorDimerMin.setBox(box);
            activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
            integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
            getController().removeAction(activityIntegrateMD);
            getController().removeAction(activityIntegrateDimer);
            getController().addAction(activityIntegrateMin);
        }

    //SET MOVABLE ATOMS
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        IAtomSet loopSet = box.getMoleculeList(fixed);
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(adAtomPos,((IAtomPositioned)((IMolecule)loopSet.getAtom(i)).getChildList().getAtom(0)).getPosition());
            if(rij.x(0)<2.0){
               movableList.add(loopSet.getAtom(i));
            } 
        }
       	for (int i=0; i<movableList.getAtomCount(); i++){
           ((IAtomPositioned)box.addNewMolecule(movable).getChildList().getAtom(0)).getPosition().E(((IAtomPositioned)((IMolecule)movableList.getAtom(i)).getChildList().getAtom(0)).getPosition());
           box.removeMolecule((IMolecule)movableList.getAtom(i));
       	}
       	movableSet = box.getMoleculeList(movable);

        //GENERATE CONFIGURATIONS FROM GAUSSIAN 
       	if(false){
	   		getController().removeAction(activityIntegrateMD);
	   		getController().removeAction(activityIntegrateDimer);
	     	RandomNumberGenerator random = new RandomNumberGenerator();
	     	IVector workVector = space.makeVector();
	   		IVector [] currentPos = new IVector [movableSet.getAtomCount()];
	   		for(int i=0; i<currentPos.length; i++){
	   			currentPos[i] = space.makeVector();
	   			currentPos[i].E(((IAtomPositioned)((IMolecule)movableSet.getAtom(i)).getChildList().getAtom(0)).getPosition());
	   		}
	     	ConfigurationFile configRandom = new ConfigurationFile(fileName+"_A_minimum");
	   		configRandom.initializeCoordinates(box);
	   		//Create multiple configurations
	   		for(int m=0; m<50; m++){
	   			WriteConfiguration genConfig = new WriteConfiguration(space);
	   			genConfig.setBox(box);
	   			genConfig.setConfName(fileName+"_config_"+m);
	   			//Displaces atom's by at most +/-0.03 in each coordinate
	   			for(int i=0; i<movableSet.getAtomCount(); i++){
	   				IVector atomPosition = ((IAtomPositioned)((IMolecule)movableSet.getAtom(i)).getChildList().getAtom(0)).getPosition();
	   				for(int j=0; j<3; j++){
	   					workVector.setX(j,0.03*random.nextGaussian());
	   				}
	   				atomPosition.Ev1Pv2(currentPos[i],workVector);
	   			}
	   			genConfig.actionPerformed();   			
	   		}
       	}
       	
    //RUN DIMER METHOD FROM GENERATED CONFIG FILE
       	if(useConfig){
            ConfigurationFile configFileFromMD = new ConfigurationFile(fileName);
        	configFileFromMD.initializeCoordinates(box);
       	}
     
    //CALCULATE VIBRATIONAL MODES
        if(calcModes){
        	getController().removeAction(activityIntegrateMD);
        	getController().removeAction(activityIntegrateDimer);
        	String file = fileName;
		    ConfigurationFile configFile = new ConfigurationFile(file);
		    configFile.initializeCoordinates(box);
		    
		    System.out.println(file+" ***Vibrational Normal Mode Analysis***");
		    System.out.println("  -Reading in system coordinates...");
		    
		    calcGradientDifferentiable = new CalcGradientDifferentiable(box, potentialMaster, movableSet, space);
		    d = new int[movableSet.getAtomCount()*3];
		    positions = new double[d.length];
		    dForces = new double[positions.length][positions.length];
		    
		    // setup position array
		    for(int i=0; i<movableSet.getAtomCount(); i++){
		        for(int j=0; j<3; j++){
		            positions[(3*i)+j] = ((IAtomPositioned)((IMolecule)movableSet.getAtom(i)).getChildList().getAtom(0)).getPosition().x(j);
		        }
		    }
		    
		    // fill dForces array
		    for(int l=0; l<d.length; l++){
		        d[l] = 1;
		        System.arraycopy(calcGradientDifferentiable.df2(d, positions), 0, dForces[l], 0, d.length);
		        System.out.println("  -Calculating force constant row "+l+"...");
		        d[l] = 0;
		    }
		    
		    calcVibrationalModes = new CalcVibrationalModes(dForces);
		    modeSigns = new int[3];
		
		    // calculate vibrational modes and frequencies
		    System.out.println("  -Calculating lambdas...");
		    lambdas = calcVibrationalModes.getLambdas();
		    System.out.println("  -Calculating frequencies...");
		    frequencies = calcVibrationalModes.getFrequencies();
		    modeSigns = calcVibrationalModes.getModeSigns();
		    
		    System.out.println("  -Writing data...");
		    // output data
		    FileWriter writer;
		    
		    //LAMBDAS
		    try { 
		        writer = new FileWriter(file+"_lambdas");
		        for(int i=0; i<lambdas.length; i++){
		            writer.write(lambdas[i]+"\n");
		        }
		        writer.close();
		    }catch(IOException e) {
		        System.err.println("Cannot open file, caught IOException: " + e.getMessage());
		        return;
		    }
		    
		    //FREQUENCIES
		    try { 
		        writer = new FileWriter(file+"_frequencies");
		        for(int i=0; i<frequencies.length; i++){
		            writer.write(frequencies[i]+"\n");
		        }
		        writer.close();
		    }catch(IOException e) {
		        System.err.println("Cannot open file, caught IOException: " + e.getMessage());
		        return;
		    }
		    
		    //MODE INFO
		    try { 
		        writer = new FileWriter(file+"_modeSigns");
		        writer.write(modeSigns[0]+" positive modes"+"\n");
		        writer.write(modeSigns[1]+" negative modes"+"\n");
		        writer.write(modeSigns[2]+" total modes"+"\n");
		        writer.close();
		    }catch(IOException e) {
		        System.err.println("Cannot open file, caught IOException: " + e.getMessage());
		        return;
		    }
		
		    System.out.println("Done.");
		    
        }
    }

    public static void main(String[] args){
    	final String APP_NAME = "DimerLJadatom";
    	final SimDimerLJadatom sim = new SimDimerLJadatom("lj_config_0", true, true, false, false, false, false);

    	sim.activityIntegrateMD.setMaxSteps(0);
    	sim.activityIntegrateDimer.setMaxSteps(1000);
    	    	
        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
        energyMeter.setBox(sim.box);
        
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorAverageCollapsing accumulatorAveragePE = new AccumulatorAverageCollapsing();
        
        DataPump energyPump = new DataPump(energyMeter,accumulatorAveragePE);       
        accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{StatType.MOST_RECENT});
        
        DisplayPlot plotPE = new DisplayPlot();
        plotPE.setLabel("PE Plot");
        
        energyAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
        accumulatorAveragePE.setPushInterval(1);      
    	
    	SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME,sim.space);
    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.add(plotPE);
    	
        sim.integratorMD.addIntervalAction(energyPump);
        sim.integratorMD.addIntervalAction(simGraphic.getPaintAction(sim.box));
        
        sim.integratorDimer.addIntervalAction(energyPump);
    	sim.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	//sim.integratorDimerMin.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
    	colorScheme.setColor(sim.fixed.getLeafType(),java.awt.Color.blue);
        colorScheme.setColor(sim.movable.getLeafType(),java.awt.Color.PINK);
    	
        simGraphic.makeAndDisplayFrame(APP_NAME);
    }

}

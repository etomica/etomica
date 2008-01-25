package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
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
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsingAverage;

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
    public Box box;
    public IVector [] saddle, normal;
    public SpeciesSpheresMono adatom, fixed, movable;
    public P2LennardJones potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    CalcGradientDifferentiable calcGradientDifferentiable;
    CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerLJadatom";
    	final SimDimerLJadatom sim = new SimDimerLJadatom("adatom", false, false, false, false);

    	sim.activityIntegrateMD.setMaxSteps(800);
    	sim.activityIntegrateDimer.setMaxSteps(700);
    	
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
        
       // MeterEnergy meterE = new MeterEnergy(sim.potentialMaster);
       // meterE.setBox(sim.box);
        
       // AccumulatorHistory totalEAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
       // DataPump totalEPump = new DataPump(meterE, totalEAccumulator);
       // totalEAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
        
    	
    	SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.add(plotPE);
    	
        sim.integratorMD.addIntervalAction(energyPump);
      //  sim.integratorMD.addIntervalAction(totalEPump);
        sim.integratorMD.addIntervalAction(simGraphic.getPaintAction(sim.box));
        
        sim.integratorDimer.addIntervalAction(energyPump);
    	sim.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
    	colorScheme.setColor(sim.fixed.getLeafType(),java.awt.Color.blue);
    	colorScheme.setColor(sim.adatom.getLeafType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getLeafType(),java.awt.Color.PINK);
    	
        simGraphic.makeAndDisplayFrame(APP_NAME);
    }

    public SimDimerLJadatom(String fileName, Boolean saddleFine, Boolean calcModes, Boolean minSearch, Boolean normalDir) {
    	super(Space3D.getInstance(), true);
    	    	
    	potentialMaster = new PotentialMaster(space);
    	
    //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
        
        BoxImposePbc imposePbc = new BoxImposePbc(box);
        
    	
    // INTEGRATOR - MD
    	integratorMD = new IntegratorVelocityVerlet(this, potentialMaster);
    	integratorMD.setTimeStep(0.01);
    	integratorMD.setTemperature(0.1);
    	integratorMD.setThermostatInterval(100);
    	integratorMD.setIsothermal(true);
    	integratorMD.setBox(box);
    	activityIntegrateMD = new ActivityIntegrate(integratorMD);
    	integratorMD.addIntervalAction(imposePbc);
    	
    //SPECIES
    	double sigma = 1.0;
    	adatom = new SpeciesSpheresMono(this);
    	Tin tinFixed = new Tin("SnFixed", Double.POSITIVE_INFINITY);
    	fixed = new SpeciesSpheresMono(this, tinFixed);
        movable = new SpeciesSpheresMono(this);      
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(adatom);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)adatom.getLeafType()).setDiameter(sigma); 
        ((AtomTypeSphere)fixed.getLeafType()).setDiameter(sigma);
        ((AtomTypeSphere)movable.getLeafType()).setDiameter(sigma);
    	
        // Must be in same order as the respective species is added to SpeciesManager
        box.setNMolecules(fixed, 256);    	
    	
    	box.setDensity(1);
    	
    	potential = new P2LennardJones(space, sigma, 1.0);
		potentialMaster.addPotential(potential, new AtomType[]{adatom.getLeafType(), adatom.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{fixed.getLeafType(), fixed.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{adatom.getLeafType(), movable.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{adatom.getLeafType(), fixed.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{movable.getLeafType(), fixed.getLeafType()});
		potentialMaster.addPotential(potential, new AtomType[]{movable.getLeafType(), movable.getLeafType()});
        
    //CRYSTAL
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc());
        config.initializeCoordinates(box); 
       
    //ADATOM CREATION AND PLACEMENT
        IMolecule iMolecule = (IMolecule)adatom.getMoleculeFactory().makeAtom();
        box.addMolecule(iMolecule);
        IVector adAtomPos = ((IAtomPositioned)iMolecule.getChildList().getAtom(0)).getPosition();
        adAtomPos.setX(0, box.getBoundary().getDimensions().x(0)/2+0.5);
        adAtomPos.setX(1, box.getBoundary().getDimensions().x(0)/16 + 0.1);
        adAtomPos.setX(2, box.getBoundary().getDimensions().x(0)/16);
        
  //INTEGRATOR - Dimer
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{adatom,movable}, fileName);
        integratorDimer.setBox(box);
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
        integratorDimer.setActivityIntegrate(activityIntegrateDimer);

    //ADD CONTROLLER ACTIONS
    	getController().addAction(activityIntegrateMD);
    	getController().addAction(activityIntegrateDimer);
    	
    //FINE-DIMER SETTINGS
        if(saddleFine==true){
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
        if(minSearch==true){
        	ConfigurationFile configFile = new ConfigurationFile(fileName+"_fine_saddle");
        	configFile.initializeCoordinates(box);
            integratorDimerMin = new IntegratorDimerMin(this, potentialMaster, new Species[]{adatom,movable}, fileName, normalDir);
            integratorDimerMin.setBox(box);
            activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
            integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
            getController().removeAction(activityIntegrateMD);
            getController().removeAction(activityIntegrateDimer);
            getController().addAction(activityIntegrateMin);
        }

        
        
    //SET MOVABLE ATOMS - top 4 layers
        
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        AtomSet loopSet = box.getMoleculeList(fixed);
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
        
    
    //CALCULATE VIBRATIONAL MODES
        if(calcModes==true){
        	getController().removeAction(activityIntegrateMD);
        	getController().removeAction(activityIntegrateDimer);
        	String file = fileName;
		    ConfigurationFile configFile = new ConfigurationFile(file);
		    configFile.initializeCoordinates(box);
		    System.out.println(file+" ***Vibrational Normal Mode Analysis***");
		    System.out.println("  -Reading in system coordinates...");
		    
		    calcGradientDifferentiable = new CalcGradientDifferentiable(box, potentialMaster, box.getLeafList().getAtomCount()-1, box.getLeafList().getAtomCount()-1);
		    d = new int[3];
		    positions = new double[d.length];
		    dForces = new double[positions.length][positions.length];
		    
		    // setup position array
		    for(int i=0; i<positions.length/3; i++){
		        for(int j=0; j<3; j++){
		            positions[(3*i)+j] = ((IAtomPositioned)box.getLeafList().getAtom(216)).getPosition().x(j);
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
}

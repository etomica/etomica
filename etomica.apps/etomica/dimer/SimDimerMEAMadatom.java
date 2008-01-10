package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
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
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.HistoryCollapsingAverage;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class SimDimerMEAMadatom extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerMEAMadatomSn";
    public final PotentialMaster potentialMaster;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorDimerMin;
    public Box box;
    public IVector [] saddle, normal;
    public SpeciesSpheresMono sn, snFix, snAdatom, ag, agFix, agAdatom, cu, cuFix, cuAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    CalcGradientDifferentiable calcGradientDifferentiable;
    CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerMEAMadatomSn";
    	final SimDimerMEAMadatom sim = new SimDimerMEAMadatom("snAdatom", false, false, false, false);

    	sim.activityIntegrateMD.setMaxSteps(7);
    	sim.activityIntegrateDimer.setMaxSteps(1);
    	sim.activityIntegrateMin.setMaxSteps(0);
    	
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
    	
    	SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.add(/*"PE Plot",*/plotPE);
    	
        sim.integratorMD.addIntervalAction(energyPump);
        sim.integratorMD.addIntervalAction(simGraphic.getPaintAction(sim.box));
        
        sim.integratorDimer.addIntervalAction(energyPump);
    	sim.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
    	
    	//Sn
    	colorScheme.setColor(sim.sn.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.snFix.getMoleculeType(),java.awt.Color.blue);
    	colorScheme.setColor(sim.snAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
        /**
        //Ag
        colorScheme.setColor(sim.ag.getMoleculeType(),java.awt.Color.darkGray);
        colorScheme.setColor(sim.agFix.getMoleculeType(),java.awt.Color.green);
        colorScheme.setColor(sim.agAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
         */
        /**
        //Cu
        colorScheme.setColor(sim.cu.getMoleculeType(),java.awt.Color.yellow);
        colorScheme.setColor(sim.cuFix.getMoleculeType(),java.awt.Color.cyan);
        colorScheme.setColor(sim.cuAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
         */
    	simGraphic.makeAndDisplayFrame(APP_NAME);
    }

    public SimDimerMEAMadatom(String fileName, Boolean saddleFine, Boolean calcModes, Boolean minSearch, Boolean normalDir) {
    	super(Space3D.getInstance(), true);
    	    	
    	potentialMaster = new PotentialMaster(space);
    	
    //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
    	
    // INTEGRATOR - MD
    	integratorMD = new IntegratorVelocityVerlet(this, potentialMaster);
    	integratorMD.setTimeStep(0.001);
    	integratorMD.setTemperature(Kelvin.UNIT.toSim(295));
    	integratorMD.setThermostatInterval(100);
    	integratorMD.setIsothermal(true);
    	integratorMD.setBox(box);
    	activityIntegrateMD = new ActivityIntegrate(integratorMD);

    //SPECIES
    	// Sn
    	Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);
    	snFix = new SpeciesSpheresMono(this, tinFixed);
        sn = new SpeciesSpheresMono(this, Tin.INSTANCE);
        snAdatom = new SpeciesSpheresMono(this, Tin.INSTANCE);
        movable = new SpeciesSpheresMono(this, Tin.INSTANCE);      
        getSpeciesManager().addSpecies(snFix);
        getSpeciesManager().addSpecies(sn);
        getSpeciesManager().addSpecies(snAdatom);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)snFix.getMoleculeType()).setDiameter(3.022); 
        ((AtomTypeSphere)sn.getMoleculeType()).setDiameter(3.022);
        ((AtomTypeSphere)snAdatom.getMoleculeType()).setDiameter(3.022);
        ((AtomTypeSphere)movable.getMoleculeType()).setDiameter(3.022);
        box.setNMolecules(snFix, 72);
    	box.setNMolecules(sn, 144);	
    	box.setNMolecules(snAdatom, 0);
    	potential = new PotentialMEAM(space);
    	potential.setParameters(snFix, ParameterSetMEAM.Sn);
		potential.setParameters(sn, ParameterSetMEAM.Sn);
		potential.setParameters(snAdatom, ParameterSetMEAM.Sn);
		potential.setParameters(movable, ParameterSetMEAM.Sn);
		this.potentialMaster.addPotential(potential, new Species[]{sn, snFix, snAdatom, movable});
        /**
        //Ag
        Silver silverFixed = new Silver("AgFix", Double.POSITIVE_INFINITY);
        agFix = new SpeciesSpheresMono(this, silverFixed);
        ag = new SpeciesSpheresMono(this, silverFixed);
        agAdatom = new SpeciesSpheresMono(this, Silver.INSTANCE);
        movable = new SpeciesSpheresMono(this, Silver.INSTANCE);
        getSpeciesManager().addSpecies(agFix);
        getSpeciesManager().addSpecies(ag);
        getSpeciesManager().addSpecies(agAdatom);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)agFix.getMoleculeType()).setDiameter(2.8895); 
        ((AtomTypeSphere)ag.getMoleculeType()).setDiameter(2.8895); 
        ((AtomTypeSphere)agAdatom.getMoleculeType()).setDiameter(2.8895);
        ((AtomTypeSphere)movable.getMoleculeType()).setDiameter(2.8895);
        box.setNMolecules(agFix, 64);
        box.setNMolecules(ag, 128);
        box.setNMolecules(agAdatom, 0);
        potential = new PotentialMEAM(space);
        potential.setParameters(agFix, ParameterSetMEAM.Ag);
        potential.setParameters(ag, ParameterSetMEAM.Ag);
        potential.setParameters(agAdatom, ParameterSetMEAM.Ag);
        potential.setParameters(movable, ParameterSetMEAM.Ag);
        this.potentialMaster.addPotential(potential, new Species[]{ag, agFix, agAdatom, movable});
         */
        /**
        //Cu
        Copper copperFixed = new Copper("CuFix", Double.POSITIVE_INFINITY);
        cuFix = new SpeciesSpheresMono(this, copperFixed);
        cu = new SpeciesSpheresMono(this, copperFixed);
        cuAdatom = new SpeciesSpheresMono(this, Copper.INSTANCE);
        movable = new SpeciesSpheresMono(this, Copper.INSTANCE);
        getSpeciesManager().addSpecies(cuFix);
        getSpeciesManager().addSpecies(cu);
        getSpeciesManager().addSpecies(cuAdatom);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)cuFix.getMoleculeType()).setDiameter(2.5561); 
        ((AtomTypeSphere)cu.getMoleculeType()).setDiameter(2.5561); 
        ((AtomTypeSphere)cuAdatom.getMoleculeType()).setDiameter(2.5561);
        ((AtomTypeSphere)movable.getMoleculeType()).setDiameter(2.5561);
        box.setNMolecules(cuFix, 64);
        box.setNMolecules(cu, 128);
        box.setNMolecules(cuAdatom, 0);
        potential = new PotentialMEAM(space);
        potential.setParameters(cuFix, ParameterSetMEAM.Cu);
        potential.setParameters(cu, ParameterSetMEAM.Cu);
        potential.setParameters(cuAdatom, ParameterSetMEAM.Cu);
        potential.setParameters(movable, ParameterSetMEAM.Cu);
        this.potentialMaster.addPotential(potential, new Species[]{cu, cuFix, cuAdatom, movable});
         */
               	
    //CRYSTAL
    	/**
    	beta-Sn box
        The dimensions of the simulation box must be proportional to those of
        the unit cell to prevent distortion of the lattice.  
        */
    	/**
    	//The values for the lattice parameters for tin's beta box (a = 5.8314 A, c = 3.1815 A) are taken from the ASM Handbook. 
    	box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
    	PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
    	BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        */
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        /**
        //Ag
        box.setDimensions(new Vector3D(4.0863*4, 4.0863*4, 4.0863*4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
         */
        /**
        //Cu
        box.setDimensions(new Vector3D(3.6148*3, 3.6148*4, 3.6148*4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
         */
        Configuration config = new ConfigurationLattice(crystal);
        config.initializeCoordinates(box); 
       
    //ADATOM CREATION AND PLACEMENT
        // Sn
        IAtom iAtom = snAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(snAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 10.0);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.1);
        ((IAtomPositioned)iAtom).getPosition().setX(2, -0.1);
        /**
        //Ag
        IAtom iAtom = agAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(agAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 7.5);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.9477016722828758);
        ((IAtomPositioned)iAtom).getPosition().setX(2, 1.0709520701043456);
        */
        /**
        //Cu
        IAtom iAtom = cuAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(cuAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 6.0);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.9477016722828758);
        ((IAtomPositioned)iAtom).getPosition().setX(2, 1.0709520701043456);
        */
        
  //INTEGRATOR - Dimer
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{snAdatom,movable}, fileName);
    	/**
    	//Ag
    	integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{agAdatom}, fileName);
    	 */
        /**
        //Cu
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{cuAdatom}, fileName);
         */
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
            integratorDimerMin = new IntegratorDimerMin(this, potentialMaster, new Species[]{snAdatom,movable}, fileName, normalDir);
            integratorDimerMin.setBox(box);
            activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
            integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
            getController().removeAction(activityIntegrateMD);
            getController().removeAction(activityIntegrateDimer);
            getController().addAction(activityIntegrateMin);
        }

        
        
    //SET MOVABLE ATOMS
        /**
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        AtomSet loopSet = box.getMoleculeList(sn);
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(((IAtomPositioned)iAtom).getPosition(),((IAtomPositioned)loopSet.getAtom(i)).getPosition()); 
            if((rij.squared())<38.0){
               movableList.add(loopSet.getAtom(i));
            } 
        }
       	for (int i=0; i<movableList.getAtomCount(); i++){
           ((IAtomPositioned)box.addNewMolecule(movable)).getPosition().E(((IAtomPositioned)movableList.getAtom(i)).getPosition());
           box.removeMolecule(movableList.getAtom(i));
       	}
        */
    
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

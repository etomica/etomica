package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;

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
import etomica.dimer.IntegratorDimerRT.PotentialMasterListDimer;
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
import etomica.nbr.CriterionSimple;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.RandomNumberGenerator;

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
    public IBox box;
    public IVector [] saddle, normal;
    public SpeciesSpheresMono fixed, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    public CalcGradientDifferentiable calcGradientDifferentiable;
    public CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    public IAtomSet movableSet;
    //public Boolean saddleFine, calcModes, minSearch, normalDir;
    

    public SimDimerMEAMadatom(String fileName, Boolean useConfig, Boolean ortho, Boolean saddleFine, Boolean calcModes, Boolean minSearch, Boolean normalDir) {
        super(Space3D.getInstance(), true);    	
    	potentialMaster = new PotentialMaster(space);
    	
    //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(random, 0, 5, space), space);
        addBox(box);
    	
    // INTEGRATOR - MD
    	integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, space);
    	integratorMD.setTimeStep(0.001);
    	integratorMD.setTemperature(Kelvin.UNIT.toSim(295));
    	integratorMD.setThermostatInterval(100);
    	integratorMD.setIsothermal(true);
    	integratorMD.setBox(box);
    	activityIntegrateMD = new ActivityIntegrate(integratorMD);

    //SPECIES
    	// Sn
    	Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);
    	fixed = new SpeciesSpheresMono(this, space, tinFixed);
        movable = new SpeciesSpheresMono(this, space, Tin.INSTANCE);      
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)fixed.getLeafType()).setDiameter(3.022); 
        ((AtomTypeSphere)movable.getLeafType()).setDiameter(3.022);
        box.setNMolecules(fixed, 420);
        
    	potential = new PotentialMEAM(space);
    	potential.setParameters(fixed.getLeafType(), ParameterSetMEAM.Sn);
		potential.setParameters(movable.getLeafType(), ParameterSetMEAM.Sn);
		
		this.potentialMaster.addPotential(potential, new AtomType[]{fixed.getLeafType(), movable.getLeafType()});
		//potentialMaster.setSpecies(new Species [] {movable});
		//potentialMaster.setRange(potential.getRange()*1.1);
		//potentialMaster.setCriterion(potential, new CriterionSimple(this, potential.getRange(), potential.getRange()*1.1));
        
        //integratorMD.addNonintervalListener(potentialMaster.getNeighborManager(box));
        //integratorMD.addIntervalAction(potentialMaster.getNeighborManager(box));  
      
		
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
        box.setDimensions(new Vector3D(5.92*3, 5.92*5, 3.23*7));
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
        Configuration config = new ConfigurationLattice(crystal, space);
        config.initializeCoordinates(box); 
       
    //ADATOM CREATION AND PLACEMENT
        // Sn
        IMolecule iMolecule = movable.makeMolecule();
        box.addMolecule(iMolecule);
        IVector adAtomPos = ((IAtomPositioned)iMolecule.getChildList().getAtom(0)).getPosition();
        adAtomPos.setX(0, 10.0);
        adAtomPos.setX(1, 0.1);
        adAtomPos.setX(2, -0.1);
        IVector newBoxLength = space.makeVector();
        newBoxLength.E(box.getBoundary().getDimensions());
        newBoxLength.setX(0, 2.0*adAtomPos.x(0)+1.0);
        box.setDimensions(newBoxLength);
        
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
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new ISpecies[]{movable}, space);
    	/**
    	//Ag
    	integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{agAdatom}, fileName);
    	 */
        /**
        //Cu
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{cuAdatom}, fileName);
         */
        //integratorDimer.addNonintervalListener(potentialMaster.getNeighborManager(box));
        //integratorDimer.addIntervalAction(potentialMaster.getNeighborManager(box));    
        integratorDimer.setBox(box);
        integratorDimer.setOrtho(ortho, false, false);
        integratorDimer.setFileName(fileName);
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
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
            if(Math.abs(rij.x(0))<5.0 && Math.abs(rij.x(1))<5.0 && Math.abs(rij.x(2))<5.0){
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

        final SimDimerMEAMadatom sim = new SimDimerMEAMadatom("meam", false, false, false, false, false, false);

        sim.activityIntegrateMD.setMaxSteps(0);
        sim.activityIntegrateDimer.setMaxSteps(0);
                
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
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, sim.space);
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.add(plotPE);
        
        sim.integratorMD.addIntervalAction(energyPump);
        sim.integratorMD.addIntervalAction(simGraphic.getPaintAction(sim.box));
        
        sim.integratorDimer.addIntervalAction(energyPump);
        sim.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
    //    sim.integratorDimerMin.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
    	
    	//Sn
    	colorScheme.setColor(sim.fixed.getLeafType(),java.awt.Color.gray);
        colorScheme.setColor(sim.movable.getLeafType(),java.awt.Color.PINK);
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

}

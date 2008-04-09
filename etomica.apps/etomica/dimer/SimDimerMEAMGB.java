package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.WriteConfiguration;
import etomica.action.XYZWriter;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.box.Box;
import etomica.chem.elements.Tin;
import etomica.config.ConfigurationFile;
import etomica.config.GrainBoundaryTiltConfiguration;
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
import etomica.space3d.Space3D;
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

public class SimDimerMEAMGB extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerMEAMadatomSn";
    public final PotentialMaster potentialMaster;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorDimerMin;
    public IBox box;
    public double clusterRadius;
    public IVector [] saddle;
    public SpeciesSpheresMono fixed, movable, dimer;
    public PotentialMEAM potential;
    public PotentialCalculationForcePressureSumGB pcGB;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    public CalcGradientDifferentiable calcGradientDifferentiable;
    public CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    public IVector adAtomPos;
    public IAtomSet movableSet;
    public int [] millerPlane;
    

    
    public SimDimerMEAMGB(String file, int[] millerPlane) {
    	super(Space3D.getInstance(), true);
        potentialMaster = new PotentialMaster(space);
        
      //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(random, 2, 5, space), space);
        addBox(box);
     
      //SPECIES
        //Sn
        Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);  
        fixed = new SpeciesSpheresMono(this, space, Tin.INSTANCE);
        movable = new SpeciesSpheresMono(this, space, Tin.INSTANCE);
        dimer = new SpeciesSpheresMono(this, space, Tin.INSTANCE);
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        getSpeciesManager().addSpecies(dimer);
        ((IAtomTypeSphere)fixed.getLeafType()).setDiameter(3.022); 
        ((IAtomTypeSphere)movable.getLeafType()).setDiameter(3.022);
        ((IAtomTypeSphere)dimer.getLeafType()).setDiameter(3.022);
        potential = new PotentialMEAM(space);
        potential.setParameters(fixed.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(movable.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(dimer.getLeafType(), ParameterSetMEAM.Sn);
        
        this.potentialMaster.addPotential(potential, new IAtomTypeLeaf[]{fixed.getLeafType(), movable.getLeafType(), dimer.getLeafType()});
                
        //Ag
        /**
        Silver silverFixed = new Silver("AgFix", Double.POSITIVE_INFINITY);
        fixed = new SpeciesSpheresMono(this, Silver.INSTANCE);
        movable = new SpeciesSpheresMono(this, Silver.INSTANCE);
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)fixed.getLeafType()).setDiameter(2.8895); 
        ((AtomTypeSphere)movable.getLeafType()).setDiameter(2.8895);
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
        potential = new PotentialMEAM(space);
        potential.setParameters(cuFix, ParameterSetMEAM.Cu);
        potential.setParameters(cu, ParameterSetMEAM.Cu);
        potential.setParameters(cuAdatom, ParameterSetMEAM.Cu);
        potential.setParameters(movable, ParameterSetMEAM.Cu);
        this.potentialMaster.addPotential(potential, new Species[]{cu, cuFix, cuAdatom, movable});
        */
    
      //CRYSTAL
    	/**
    	Sn
    	beta-Sn box
        
        The dimensions of the simulation box must be proportional to those of
        the unit cell to prevent distortion of the lattice.  The values for the 
        lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        angstroms) are taken from the ASM Handbook. 
    	*/
    	double a = 5.92; 
    	double c = 3.23;
    	//box.setDimensions(new Vector3D((Math.sqrt( (4*Math.pow(c, 2))+Math.pow(a,2)))*3, a*3, c*10));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, a, c);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new ISpecies[] {fixed, movable}, 4.5, space);
        /**
        //Ag
        box.setDimensions(new Vector3D(4.0863*4, 4.0863*4, 4.0863*4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new Species[] {fixed, movable}, 4.56, space);
        */
        /**
        //Cu
        box.setDimensions(new Vector3D(3.6148*5, 3.6148*4, 3.6148*4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new Species[] {cuFix, cu}, 4.56);
        */

        gbtilt.setFixedSpecies(fixed);
        gbtilt.setMobileSpecies(movable);
        gbtilt.setGBplane(millerPlane);
        gbtilt.setBoxSize(box, new int[] {4,4,10});
        gbtilt.initializeCoordinates(box);
    }
    
    public void setMovableAtoms(double distance, IVector center){
        distance = distance*distance;
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        IAtomSet loopSet = box.getMoleculeList(fixed);
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(center,((IAtomPositioned)((IMolecule)loopSet.getAtom(i)).getChildList().getAtom(0)).getPosition());
            if(rij.squared() < distance){
               movableList.add(loopSet.getAtom(i));
            } 
        }
        for (int i=0; i<movableList.getAtomCount(); i++){
           ((IAtomPositioned)box.addNewMolecule(movable).getChildList().getAtom(0)).getPosition().E(((IAtomPositioned)((IMolecule)movableList.getAtom(i)).getChildList().getAtom(0)).getPosition());
           box.removeMolecule((IMolecule)movableList.getAtom(i));
        }
        movableSet = box.getMoleculeList(movable);
    }
    
    public void initializeConfiguration(String fileName){
        ConfigurationFile config = new ConfigurationFile(fileName);
        config.initializeCoordinates(box);
    }
    
    public void generateConfigs(String fileName, double percentd){       
        
        RandomNumberGenerator random = new RandomNumberGenerator();
        IVector workVector = space.makeVector();
        IVector [] currentPos = new IVector [movableSet.getAtomCount()];
        for(int i=0; i<currentPos.length; i++){
            currentPos[i] = space.makeVector();
            currentPos[i].E(((IAtomPositioned)((IMolecule)movableSet.getAtom(i)).getChildList().getAtom(0)).getPosition());
        }
        
        //Create multiple configurations
        for(int m=0; m<50; m++){
            WriteConfiguration genConfig = new WriteConfiguration(space);
            genConfig.setBox(box);
            genConfig.setConfName(fileName+"_config_"+m);
            //Displaces atom's by at most +/-0.03 in each coordinate
            for(int i=0; i<movableSet.getAtomCount(); i++){
                IVector atomPosition = ((IAtomPositioned)((IMolecule)movableSet.getAtom(i)).getChildList().getAtom(0)).getPosition();
                for(int j=0; j<3; j++){
                    workVector.setX(j,percentd*random.nextGaussian());
                }
                atomPosition.Ev1Pv2(currentPos[i],workVector);
            }
            genConfig.actionPerformed();            
        }
    }
    
    public void calculateVibrationalModes(String fileName){
        
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
    
    public void enableMolecularDynamics(long maxSteps){
        integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integratorMD.setTimeStep(0.001);
        integratorMD.setTemperature(Kelvin.UNIT.toSim(0));
        integratorMD.setThermostatInterval(100);
        integratorMD.setIsothermal(true);
        integratorMD.setBox(box);
        pcGB = new PotentialCalculationForcePressureSumGB(space, box);
        integratorMD.setForceSum(pcGB);
        //integratorMD.addNonintervalListener(potentialMaster.getNeighborManager(box));
       // integratorMD.addIntervalAction(potentialMaster.getNeighborManager(box));  
        activityIntegrateMD = new ActivityIntegrate(integratorMD);
        getController().addAction(activityIntegrateMD);
        activityIntegrateMD.setMaxSteps(maxSteps);
    }
    
    public void enableDimerSearch(String fileName, long maxSteps, Boolean orthoSearch, Boolean fine){
        
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new ISpecies[]{movable}, space);
        //integratorDimer.addNonintervalListener(potentialMaster.getNeighborManager(box));
        //integratorDimer.addIntervalAction(potentialMaster.getNeighborManager(box));    
        integratorDimer.setBox(box);
        integratorDimer.setOrtho(orthoSearch, false);
        if(fine){
            ConfigurationFile configFile = new ConfigurationFile(fileName+"_saddle");
            configFile.initializeCoordinates(box);
            
            integratorDimer.setFileName(fileName+"_fine");
            integratorDimer.deltaR = 0.0005;
            integratorDimer.dXl = 10E-5;       
            integratorDimer.deltaXmax = 0.005;
            integratorDimer.dFsq = 0.0001*0.0001;
            integratorDimer.dFrot = 0.01;
        }
        integratorDimer.setFileName(fileName);
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
        integratorDimer.setActivityIntegrate(activityIntegrateDimer);
        getController().addAction(activityIntegrateDimer);
        activityIntegrateDimer.setMaxSteps(maxSteps);
    }
        
    public void enableMinimumSearch(String fileName, Boolean normalDir, Boolean onlyMin){
        
        if(onlyMin){
            ConfigurationFile configFile = new ConfigurationFile(fileName);
            configFile.initializeCoordinates(box);
        }
        else{
            fileName = fileName+"_fine";
            ConfigurationFile configFile = new ConfigurationFile(fileName);
            configFile.initializeCoordinates(box);
        }
        integratorDimerMin = new IntegratorDimerMin(this, potentialMaster, new ISpecies[]{movable}, fileName, normalDir, space);
        integratorDimerMin.setBox(box);
        activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
        integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
        getController().addAction(activityIntegrateMin);
    }
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerMEAMadatomGB";
    	final SimDimerMEAMGB sim = new SimDimerMEAMGB("sngb", new int[] {1,1,1});
    	
    	sim.activityIntegrateMD.setMaxSteps(900);
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
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1, sim.space);
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName("sngb");
        xyzwriter.setIsAppend(true);
        
        simGraphic.add(/*"PE Plot",*/plotPE);
        /*
        //Load in MD minimized configuration
        ConfigurationFile configurationFile = new ConfigurationFile("filename");
        configurationFile.initializeCoordinates(sim.box);
        
        sim.integratorMD.addIntervalAction(energyPump);
        sim.integratorMD.addIntervalAction(simGraphic.getPaintAction(sim.box));
        sim.integratorMD.addIntervalAction(xyzwriter);
        sim.integratorMD.setActionInterval(xyzwriter, 10);
        
        sim.integratorDimer.addIntervalAction(energyPump);
        sim.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 10);
        */
        
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        
        //Sn
        colorScheme.setColor(sim.fixed.getLeafType(),java.awt.Color.blue);
        colorScheme.setColor(sim.movable.getLeafType(),java.awt.Color.gray);
        colorScheme.setColor(sim.dimer.getLeafType(), java.awt.Color.white);
        
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

package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
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
import etomica.units.ElectronVolt;
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
    public Box box;
    public IVector [] saddle;
    public SpeciesSpheresMono sn, snFix, snAdatom, ag, agFix, agAdatom, cu, cuFix, cuAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer;
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerMEAMadatomSn";
    	final SimDimerMEAMadatom sim = new SimDimerMEAMadatom();
    	
    	sim.activityIntegrateMD.setMaxSteps(7);
    	sim.activityIntegrateDimer.setMaxSteps(1);
    	
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

    public SimDimerMEAMadatom() {
    	super(Space3D.getInstance(), true);
    	
    	potentialMaster = new PotentialMaster(space);
    	
    	integratorMD = new IntegratorVelocityVerlet(this, potentialMaster);
    	
    	integratorMD.setTimeStep(0.001);
    	integratorMD.setTemperature(Kelvin.UNIT.toSim(295));
    	integratorMD.setThermostatInterval(100);
    	integratorMD.setIsothermal(true);
    	
    	activityIntegrateMD = new ActivityIntegrate(integratorMD);
    	
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
         */
        
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
        
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{snAdatom}, "SnAdatom");
    	/**
    	//Ag
    	integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{agAdatom}, "AgAdatom");
    	 */
        /**
        //Cu
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{cuAdatom}, "CuAdatom");
         */
        
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);

    	// First simulation style
    	getController().addAction(activityIntegrateMD);
    	// Second simulation style
    	getController().addAction(activityIntegrateDimer);
    	
    	// Sn
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
    	
    	integratorMD.setBox(box);
    	integratorDimer.setBox(box);	

    	// beta-Sn box
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        
    	//box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
    	//PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
    	//BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        
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
       
       
       
       
    }
}

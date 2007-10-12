package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.GrainBoundaryTiltConfiguration;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
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

public class SimDimerMEAMadatomCuGB extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM Md3D";
    public final PotentialMaster potentialMaster;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorMEP;
    public Box box;
    public IVector [] saddle;
    public SpeciesSpheresMono sn, snFix, snAdatom, cu, cuFix, cuAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMEP;
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerMEAMadatomCu";
    	final SimDimerMEAMadatomCuGB sim = new SimDimerMEAMadatomCuGB();
    	
    	sim.activityIntegrateMD.setMaxSteps(500);
    	sim.activityIntegrateDimer.setMaxSteps(1000);
    	sim.activityIntegrateMEP.setMaxSteps(250);
    	
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
        
        DataLogger elog = new DataLogger();
        elog.setFileName("cu_sim-energy2");
        elog.setAppending(true);
        elog.setWriteInterval(1);
        elog.setCloseFileEachTime(true);
        elog.setDataSink(new DataTableWriter());
        sim.getController().getEventManager().addListener(elog);
        
        DataPump pump = new DataPump(energyMeter, elog);
        
        sim.integratorDimer.addIntervalAction(pump);
        sim.integratorDimer.setActionInterval(pump, 1);
        sim.integratorDimer.addIntervalAction(energyPump);
    	sim.integratorDimer.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	
    	sim.integratorMEP.addIntervalAction(pump);
    	sim.integratorMEP.setActionInterval(pump,1);
    	sim.integratorMEP.addIntervalAction(energyPump);
    	sim.integratorMEP.addIntervalAction(simGraphic.getPaintAction(sim.box));
    	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());

    	//Cu
    	///**
    	colorScheme.setColor(sim.cu.getMoleculeType(),java.awt.Color.yellow);
        colorScheme.setColor(sim.cuFix.getMoleculeType(),java.awt.Color.yellow);
        colorScheme.setColor(sim.cuAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
        //**/
    	
    	simGraphic.makeAndDisplayFrame(APP_NAME);
    	
    }
    
    
    public SimDimerMEAMadatomCuGB() {
    	super(Space3D.getInstance(), true);
    	
    	potentialMaster = new PotentialMaster(space);
    	
    	integratorMD = new IntegratorVelocityVerlet(this, potentialMaster);
    	
    	integratorMD.setTimeStep(0.001);
    	integratorMD.setTemperature(Kelvin.UNIT.toSim(100));
    	integratorMD.setThermostatInterval(100);
    	integratorMD.setIsothermal(true);
    	
    	activityIntegrateMD = new ActivityIntegrate(integratorMD);

        // Cu
        ///**
        Copper copperFixed = new Copper("CuFix", Double.POSITIVE_INFINITY);
        
        cuFix = new SpeciesSpheresMono(this, copperFixed);
        cu = new SpeciesSpheresMono(this, Copper.INSTANCE);
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
        //*/ 
        
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
        
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{cuAdatom}, "CuGB");     
    	activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
    	        
    	integratorMEP = new IntegratorDimerMin(this, potentialMaster, new Species[]{cuAdatom, cu});
    	activityIntegrateMEP = new ActivityIntegrate(integratorMEP);
    	    	
    	// First simulation style
    	getController().addAction(activityIntegrateMD);
    	// Second simulation style
    	getController().addAction(activityIntegrateDimer);
    	// Third simulation style
    	getController().addAction(activityIntegrateMEP);
    	
        
        potential = new PotentialMEAM(space);
        
    	potential.setParameters(cuFix, ParameterSetMEAM.Cu);
        potential.setParameters(cu, ParameterSetMEAM.Cu);
        potential.setParameters(cuAdatom, ParameterSetMEAM.Cu);
        potential.setParameters(movable, ParameterSetMEAM.Cu);
        
        this.potentialMaster.addPotential(potential, new Species[]{cu, cuFix, cuAdatom, movable});
		//*/		
    	
    	integratorMD.setBox(box);
    	integratorDimer.setBox(box);	
    	integratorMEP.setBox(box);

        // Cu box
        ///**
        box.setDimensions(new Vector3D(3.6148*5, 3.6148*4, 3.6148*4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new Species[] {cuFix, cu}, 4.56);
        gbtilt.setRotation(1, 53.1*Math.PI/180);
        gbtilt.initializeCoordinates(box); 

        // Cu
        /**
        IAtom iAtom = cuAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(cuAdatom).addChildAtom(iAtom);
       
        // ((IAtomPositioned)iAtom).getPosition().setX(0, 7.912128210706072);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 6.246340066365732);
        ((IAtomPositioned)iAtom).getPosition().setX(1, 0.9477016722828758);
        ((IAtomPositioned)iAtom).getPosition().setX(2, 1.0709520701043456);
        
        // ((IAtomPositioned)iAtom).getPosition().setX(0, 5.892029514644299);
        // ((IAtomPositioned)iAtom).getPosition().setX(1, 1.8461220902169309);
        // ((IAtomPositioned)iAtom).getPosition().setX(2, 1.7526031843740997);
        
       /** 
        IVector rij = space.makeVector();
        AtomArrayList movableList = new AtomArrayList();
        AtomSet loopSet = box.getMoleculeList(cu);
        
        for (int i=0; i<loopSet.getAtomCount(); i++){
            rij.Ev1Mv2(((IAtomPositioned)iAtom).getPosition(),((IAtomPositioned)loopSet.getAtom(i)).getPosition());
          
            if(rij.squared()<21.0){
               movableList.add(loopSet.getAtom(i));
            } 
        }
       
       for (int i=0; i<movableList.getAtomCount(); i++){
           ((IAtomPositioned)box.addNewMolecule(movable)).getPosition().E(((IAtomPositioned)movableList.getAtom(i)).getPosition());
           box.removeMolecule(movableList.getAtom(i));
       }
        
        
       */
    }
}

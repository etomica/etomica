package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Tin;
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

public class DimerMEAMadatomSn extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerMEAMadatomSn";
    public final PotentialMaster potentialMaster;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMEP integratorMEP;
    public Box box;
    public IVector [] saddle;
    public SpeciesSpheresMono sn, snFix, snAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMEP;
    
    public static void main(String[] args){
    	final String APP_NAME = "DimerMEAMadatom";
    	final DimerMEAMadatomSn sim = new DimerMEAMadatomSn();
    	
    	sim.activityIntegrateMD.setMaxSteps(500);
    	sim.activityIntegrateDimer.setMaxSteps(500);
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
    	
    	//Sn

    	colorScheme.setColor(sim.sn.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.snFix.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.snAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
    	
    	simGraphic.makeAndDisplayFrame(APP_NAME);
    	
    }
    
    
    public DimerMEAMadatomSn() {
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
        
        box = new Box(new BoundaryRectangularSlit(space, random, 0, 5));
        addBox(box);
        
        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new Species[]{snAdatom, sn});     
    	activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
    	        
    	integratorMEP = new IntegratorDimerMEP(this, potentialMaster, saddle, 10E-4, 10E-4, new Species[]{snAdatom, sn});
    	activityIntegrateMEP = new ActivityIntegrate(integratorMEP);
    	    	
    	// First simulation style
    	getController().addAction(activityIntegrateMD);
    	// Second simulation style
    	getController().addAction(activityIntegrateDimer);
    	// Third simulation style
    	getController().addAction(activityIntegrateMEP);
    	
    	// Sn
    	box.setNMolecules(snFix, 64);
    	box.setNMolecules(sn, 128);	
    	box.setNMolecules(snAdatom, 0);
    	
    	potential = new PotentialMEAM(space);
		
    	potential.setParameters(snFix, ParameterSetMEAM.Sn);
		potential.setParameters(sn, ParameterSetMEAM.Sn);
		potential.setParameters(snAdatom, ParameterSetMEAM.Sn);
		potential.setParameters(movable, ParameterSetMEAM.Sn);
		
		this.potentialMaster.addPotential(potential, new Species[]{sn, snFix, snAdatom, movable});
    	
    	integratorMD.setBox(box);
    	integratorDimer.setBox(box);	
    	integratorMEP.setBox(box);
    	
    	// beta-Sn box
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        
    	box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        
        // Sn
       
        IAtom iAtom = snAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(snAdatom).addChildAtom(iAtom);
        ((IAtomPositioned)iAtom).getPosition().setX(0, 1.62*5.8314);
       
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

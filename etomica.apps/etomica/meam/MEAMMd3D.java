package etomica.meam;
import java.util.ArrayList;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.IDataInfo;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.PotentialMasterList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.CompoundUnit;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Unit;
import etomica.util.HistoryCollapsingAverage;

/**
 * Molecular-Dynamics Simulation Using the Modified Embedded-Atom Method 
 * (MEAM) Potential.  
 * 
 * The MEAM potential is intended for use with metallic and covalently-bonded
 * solid systems.
 * 
 * The MEAM potential for an atom is built using terms describing parts of the
 * relationships between the atom and each of its neighbors, the number of which 
 * is determined by a cutoff and/or screening function.  Each type of pair-
 * wise term is summed over all the neighbors, and then used in expressions 
 * describing the embedding energy and the repulsive energy of the atom.   
 * Effectively, the MEAM potential is a many-body potential.  
 * 
 * This class was adapted from LjMd3D.java by K.R. Schadel and A. Schultz in July 
 * 2005.  Intitially, it employed a version of the embedded-atom method potential, 
 * and was later adapted in February 2006 to use the modified embedded-atom method
 * potential.
 */
 
public class MEAMMd3D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM Md3D";
    public final PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono sn;
    public SpeciesSpheresMono ag;
    public SpeciesSpheresMono cu;
    public Box box;
    public PotentialMEAM potentialN;
    public Controller controller;
    public DisplayBox display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public IDataInfo info2;

    public static void main(String[] args) {
    	MEAMMd3D sim = new MEAMMd3D();
    	
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
    	MeterKineticEnergy kineticMeter = new MeterKineticEnergy();
    	
    	energyMeter.setBox(sim.box);
    	kineticMeter.setBox(sim.box);
        
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory kineticAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        
        AccumulatorAverage accumulatorAveragePE = new AccumulatorAverage(50);
    	AccumulatorAverage accumulatorAverageKE = new AccumulatorAverage(50);
    	
    	DataPump energyPump = new DataPump(energyMeter,accumulatorAveragePE);   	
    	DataPump kineticPump = new DataPump(kineticMeter, accumulatorAverageKE);
    	
    	accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{StatType.MOST_RECENT});
    	accumulatorAverageKE.addDataSink(kineticAccumulator, new StatType[]{StatType.MOST_RECENT});
    	
    	DisplayPlot plotPE = new DisplayPlot();
        plotPE.setLabel("PE Plot");
        DisplayPlot plotKE = new DisplayPlot();
        plotKE.setLabel("KE Plot");
    	
        energyAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
        kineticAccumulator.setDataSink(plotKE.getDataSet().makeDataSink());
        //energyAccumulator.setBlockSize(50);
        
    	accumulatorAveragePE.setPushInterval(1);
    	accumulatorAverageKE.setPushInterval(1);

    	//Heat Capacity (PE)
    	DataProcessorCvMD dataProcessorPE = new DataProcessorCvMD();
    	dataProcessorPE.setIntegrator(sim.integrator);
    	
    	//Heat Capacity (KE)
    	DataProcessorCvMD dataProcessorKE = new DataProcessorCvMD();
    	dataProcessorKE.setIntegrator(sim.integrator);
    	
    	accumulatorAveragePE.addDataSink(dataProcessorPE, new StatType[]{StatType.STANDARD_DEVIATION});
    	accumulatorAverageKE.addDataSink(dataProcessorKE, new StatType[]{StatType.STANDARD_DEVIATION});

        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.addIntervalAction(kineticPump);
       
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        ArrayList dataStreamPumps = simGraphic.getController().getDataStreamPumps();
        dataStreamPumps.add(energyPump);
        dataStreamPumps.add(kineticPump);
        
    	DisplayTextBox cvBoxPE = new DisplayTextBox();
    	dataProcessorPE.setDataSink(cvBoxPE);
    	cvBoxPE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
    	cvBoxPE.setLabel("PE Cv contrib.");
    	DisplayTextBox cvBoxKE = new DisplayTextBox();
    	dataProcessorKE.setDataSink(cvBoxKE);
    	cvBoxKE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
    	cvBoxKE.setLabel("KE Cv contrib.");

    	simGraphic.add(/*"PE Plot",*/plotPE);
    	simGraphic.add(/*"KE Plot",*/plotKE);
    	
    	simGraphic.getPanel().controlPanel.add(cvBoxKE.graphic(), SimulationPanel.getVertGBC());
    	simGraphic.getPanel().controlPanel.add(cvBoxPE.graphic(), SimulationPanel.getVertGBC());

    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
    	colorScheme.setColor(sim.sn.getMoleculeType(),java.awt.Color.blue);
    	colorScheme.setColor(sim.ag.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.cu.getMoleculeType(),java.awt.Color.orange);

    	simGraphic.makeAndDisplayFrame(APP_NAME);

    	//sim.activityIntegrate.setMaxSteps(1000);
    	//sim.getController().run();
    	//DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
        //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
                    // /sim.species.getAgent(sim.box).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.box).getNMolecules();  // orig line
        //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
    }
    
    public MEAMMd3D() {
        super(Space3D.getInstance(), true); //INSTANCE); kmb change 8/3/05
        potentialMaster = new PotentialMasterList(this);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(295));
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        sn = new SpeciesSpheresMono(this, Tin.INSTANCE);
        ag = new SpeciesSpheresMono(this, Silver.INSTANCE);
        cu = new SpeciesSpheresMono(this, Copper.INSTANCE);

        getSpeciesManager().addSpecies(sn);
        getSpeciesManager().addSpecies(ag);
        getSpeciesManager().addSpecies(cu);

        /** The following values come from either the ASM Handbook or Cullity & Stock's 
         * "Elements of X-Ray Diffraction" (2001)
         */
        ((AtomTypeSphere)sn.getMoleculeType()).setDiameter(3.022); 
        
        ((AtomTypeSphere)ag.getMoleculeType()).setDiameter(2.8895); 
        
        ((AtomTypeSphere)cu.getMoleculeType()).setDiameter(2.5561); 
        
        
        box = new Box(this);
        addBox(box);
        box.setNMolecules(sn, 0);
        box.setNMolecules(ag, 256);
        box.setNMolecules(cu, 0);
        
        // beta-Sn box
        /**
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        //box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        LatticeCrystal crystal = new LatticeCrystal(new Crystal(
        		primitive, new BasisBetaSnA5(primitive)));
		*/
		
        //FCC Cu
		/**
	    box.setDimensions(new Vector3D(3.6148*4, 3.6148*4, 3.6148*4));
	    PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
	    LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		        primitive, new BasisCubicFcc(primitive)));
	    */
        
        //FCC Ag
        
	    box.setDimensions(new Vector3D(4.0863*4, 4.0863*4, 4.0863*4));
	    PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
	    BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
	    
           
		Configuration config = new ConfigurationLattice(crystal);
		config.initializeCoordinates(box);
        
		potentialN = new PotentialMEAM(space);
		potentialN.setParameters(sn, ParameterSetMEAM.Sn);
		potentialN.setParameters(ag, ParameterSetMEAM.Ag);
		potentialN.setParameters(cu, ParameterSetMEAM.Cu);
		potentialN.setParametersIMC(cu, ParameterSetMEAM.Cu3Sn);
		potentialN.setParametersIMC(ag, ParameterSetMEAM.Ag3Sn);
        this.potentialMaster.addPotential(potentialN, new Species[]{sn, ag, cu});    
        potentialMaster.setRange(potentialN.getRange()*1.1);
        potentialMaster.setCriterion(potentialN, new CriterionSimple(this, potentialN.getRange(), potentialN.getRange()*1.1));
        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));
        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
        
        integrator.setBox(box);
		
        // IntegratorCoordConfigWriter - Displacement output (3/1/06 - MS)
        //IntegratorCoordConfigWriter coordWriter = new IntegratorCoordConfigWriter(space, "MEAMoutput");
        //coordWriter.setBox(box);
        //coordWriter.setIntegrator(integrator);
        //coordWriter.setWriteInterval(100);
        
        // Control simulation lengths
        //activityIntegrate.setMaxSteps(500);

		energy = new MeterEnergy(potentialMaster);
    }
    
}
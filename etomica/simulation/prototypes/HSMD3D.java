package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.ParameterBase;

/**
 * 
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.  
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class HSMD3D extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    private static final long serialVersionUID = 1L;
    /**
     * The Box holding the atoms. 
     */
    public final IBox box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2SquareWell potential;
    
    public final IPotentialMaster potentialMaster;
    
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public HSMD3D(ISpace _space) {
        this(_space, new HSMD3DParam());
    }
    
    public HSMD3D(ISpace _space, HSMD3DParam params) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(_space);

        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, 3.0, space) : new PotentialMasterMonatomic(this);

        int numAtoms = params.nAtoms;
        double neighborRangeFac = 1.2;
        double sigma = 1.0;
        double lambda = 2.5;
        if (params.useNeighborLists) {
            ((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*sigma*lambda);
        }

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setTemperature(0.18);
        integrator.setIsothermal(true);
        integrator.setTimeStep(0.02);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        potential = new P2SquareWell(space, sigma, 2.5, -1.0, false);
        IAtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new IAtomType[]{leafType, leafType});

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.eta * 2 * space.D() / Math.PI);
        inflater.setTargetDensity(0.38);
        inflater.actionPerformed();
        if (space.D() == 3) {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        }
        else {
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        }
        //deformed
//        box.setBoundary(
//            new etomica.space.BoundaryDeformablePeriodic(
//            space,getRandom(),
//            new IVector[]{
//              new Vector3D(-4,1,1),
//              new Vector3D(2,6,4),
//              new Vector3D(1,2,6)}));
        //truncated octahedron
//        box.setBoundary(
//            new etomica.space3d.BoundaryTruncatedOctahedron(this));
        
        integrator.setBox(box);

        if (params.useNeighborLists) { 
            NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager(box);
            integrator.getEventManager().addListener(nbrManager);
        }
        else {
            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        }
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
    	final String APP_NAME = "HSMD3D";

    	ISpace sp = Space3D.getInstance();
        HSMD3DParam params = new HSMD3DParam();
        final HSMD3D sim = new HSMD3D(sp, params);
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, sim.space, sim.getController());
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim, sp, sim.getController()));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);

        nSelector.setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(nSelector);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);
        
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorHistory peHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
        peHistory.setTimeDataSource(new DataSourceCountTime(sim.integrator));
        DataPumpListener pePump = new DataPumpListener(meterPE, peHistory, 100);
        sim.integrator.getEventManager().addListener(pePump);
        DisplayPlot pePlot = new DisplayPlot();
        peHistory.setDataSink(pePlot.getDataSet().makeDataSink());

        pePlot.setLabel("PE");
        simGraphic.add(pePlot);
    }

    public static HSMD3DParam getParameters() {
        return new HSMD3DParam();
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class HSMD3DParam extends ParameterBase {
        public int nAtoms = 1024;
        public double eta = 0.35;
        public boolean useNeighborLists = true;
    }
}

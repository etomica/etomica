package etomica.modules.render;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardBondedList;
import etomica.potential.P2Ideal;
import etomica.potential.P2PenetrableSquareWell;
import etomica.potential.PotentialHard;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
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
public class RenderMD extends Simulation {

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

    public final P2HardBondedList potential;
    public final PotentialHard potentialBonded, potentialNonBonded;
    
    public final IPotentialMaster potentialMaster;
    
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public RenderMD(ISpace _space) {
        this(_space, new RenderMDParam());
    }
    
    public RenderMD(ISpace _space, RenderMDParam params) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(_space);

        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, 3.0, space) : new PotentialMasterMonatomic(this);
        
        ParseObj parser = new ParseObj(params.file);

        int numAtoms = parser.nAtoms;
        
        
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
        
        potentialBonded = new P2PenetrableSquareWell(space);
        potentialNonBonded = new P2Ideal(space);
        potential = new P2HardBondedList(this, potentialBonded, potentialNonBonded);
        
        IAtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new IAtomType[]{leafType, leafType});

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.eta * 2 * space.D() / Math.PI);
        inflater.setTargetDensity(15.0);
        inflater.actionPerformed();
        
        IAtomList leafList = box.getLeafList();
        for (int iLeaf=0; iLeaf<numAtoms; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            IVectorMutable pos = a.getPosition();
            pos.E(parser.positions[iLeaf]);
        }

        
//        new ConfigurationFile(params.file).initializeCoordinates(box);
//        new ConfigurationCar(space).initializeCoordinates(box);
        
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
    	final String APP_NAME = "RenderMD";

    	ISpace sp = Space3D.getInstance();
        RenderMDParam params = new RenderMDParam();
        final RenderMD sim = new RenderMD(sp, params);
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, sim.space, sim.getController());
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim, sp, sim.getController()));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);
        ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getAtomType(0),0.1);

        nSelector.setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(nSelector);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);
    }

    public static RenderMDParam getParameters() {
        return new RenderMDParam();
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class RenderMDParam extends ParameterBase {
        public int nAtoms = 6750;
        public double eta = 0.7;
        public boolean useNeighborLists = false;
        public String file = "/Users/kofke/Documents/workspace/car self-assembly/mustang";
    }
}

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeLeaf;
import etomica.chem.models.ModelChain;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.graphics.BondListener;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.box.Box;
import etomica.potential.P2HardBond;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;

public class ChainHSMD3D extends Simulation {

    private static final long serialVersionUID = 2L;
    public Box box;
    public IntegratorHard integrator;
    public SpeciesSpheres species;
    public P2HardSphere potential;
    private ModelChain model;
    
    public ChainHSMD3D() {
        this(Space3D.getInstance());
    }
    private ChainHSMD3D(Space space) {
        super(space, true);
        PotentialMasterList potentialMaster = new PotentialMasterList(this);
        int numAtoms = 704;
        int chainLength = 4;
        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac);

        integrator = new IntegratorHard(this, potentialMaster);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator, false, true);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        
        model = new ModelChain();
        model.setNumAtoms(chainLength);
        model.setBondingPotential(new P2HardBond(space, 1.0, 0.15, true));
        
        species = (SpeciesSpheres)model.makeSpecies(this);
        potentialMaster.addModel(model);
        ((ConformationLinear)model.getConformation()).setBondLength(1.0);
        ((ConformationLinear)model.getConformation()).setAngle(1,0.5);
        
        box = new Box(this);
        double l = 14.4573*Math.pow((chainLength*numAtoms/2020.0),1.0/3.0);
        box.getBoundary().setDimensions(Space.makeVector(new double[]{l,l,l}));
        addBox(box);
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc());
        box.setNMolecules(species, numAtoms);
        config.initializeCoordinates(box);
        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));

        BoxImposePbc pbc = new BoxImposePbc(box);
        integrator.addIntervalAction(pbc);
        pbc.setApplyToMolecules(true);
        
        potential = new P2HardSphere(space, 1.0, true);
        AtomTypeLeaf leafType = (AtomTypeLeaf)((AtomTypeGroup)species.getMoleculeType()).getChildTypes()[0];
        potentialMaster.addPotential(potential, new AtomType[]{leafType,leafType});

        integrator.setBox(box);
    }

    public static void main(String[] args) {
      final etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
      final SimulationGraphic simGraphic = new SimulationGraphic(sim);
      BondListener bl = new BondListener(sim.box,(DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas);
      bl.addModel(sim.model);

      simGraphic.getController().getReinitButton().setPostAction(simGraphic.getDisplayBoxPaintAction(sim.box));

      simGraphic.makeAndDisplayFrame();
    }
}//end of class

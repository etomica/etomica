package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.atom.AtomTypeMolecule;
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
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;

public class ChainHSMD3D extends Simulation {

    private static final long serialVersionUID = 2L;
    public IBox box;
    public IntegratorHard integrator;
    public SpeciesSpheres species;
    public P2HardSphere potential;
    private ModelChain model;

    private ChainHSMD3D() {
        super(Space3D.getInstance(), true);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        int numAtoms = 108;
        int chainLength = 4;
        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac);

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator, 1, true);
        getController().addAction(activityIntegrate);

        model = new ModelChain();
        model.setNumAtoms(chainLength);
        model.setBondingPotential(new P2HardBond(space, 1.0, 0.15, true));

        species = (SpeciesSpheres)model.makeSpecies(this);
        potentialMaster.addModel(model);
        ((ConformationLinear)model.getConformation()).setBondLength(1.0);
        ((ConformationLinear)model.getConformation()).setAngle(1,0.35);
        
        box = new Box(this, space);
        double l = 14.4573*Math.pow((chainLength*numAtoms/2020.0),1.0/3.0);
        box.getBoundary().setDimensions(Space.makeVector(new double[]{l,l,l}));
        addBox(box);
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        box.setNMolecules(species, numAtoms);
        config.initializeCoordinates(box);
        integrator.addIntervalAction(potentialMaster.getNeighborManager(box));
        integrator.addNonintervalListener(potentialMaster.getNeighborManager(box));

        potential = new P2HardSphere(space, 1.0, true);
        AtomTypeLeaf leafType = species.getLeafType();
        potentialMaster.addPotential(potential, new IAtomType[]{leafType,leafType});

        integrator.setBox(box);
    }

    public static void main(String[] args) {

      final etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
      final SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space);
      BondListener bl = new BondListener(sim.box,(DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas);
      bl.addModel(sim.model);

      simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

      simGraphic.makeAndDisplayFrame();
    }
}//end of class

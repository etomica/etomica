package etomica.osmoticvirial;

import etomica.action.ActionIntegrate;
import etomica.action.BoxInflate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;

/**
 * Created by aksharag on 6/16/17.
 */
public class fep extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species1; //solvent
    public SpeciesSpheresMono species2; //solute
    public Box box;
    public P2LennardJones potential;
    public Controller controller;

    public fep(int numAtoms, int numSteps){
        super(Space3D.getInstance());
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this,space);
        double sigma = 1.0;
        integrator = new IntegratorMC(this, potentialMaster);
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        mcMoveAtom.setStepSize(0.2*sigma);
        ((MCMoveStepTracker)mcMoveAtom.getTracker()).setTunable(false);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().setEquilibrating(false);
        ActionIntegrate actionIntegrate = new ActionIntegrate(integrator, false);
        actionIntegrate.setMaxSteps(numSteps);
        getController().addAction(actionIntegrate);
        species1 = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        addSpecies(species2);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species1,numAtoms); //TODO: add for species 2
        BoxInflate inflater = new BoxInflate(box,space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1);
        double truncationRadius = 3.0*sigma;
        if(truncationRadius>0.5*box.getBoundary().getBoxSize().getX(0)){
            throw new RuntimeException(" Truncation radius is too large. Max allowed is:"+ 0.5*box.getBoundary().getBoxSize().getX(0));
        }

        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        potentialMaster.setCellRange(3);
        potentialMaster.setRange(potentialTruncated.getRange());
        AtomType leafType = species1.getLeafType(); //TODO: add for species2
        potentialMaster.addPotential(potentialTruncated, new AtomType[]{leafType, leafType});
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);

        integrator.setBox(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();


    }

    public static void main(String[] args){
        simParams params = new simParams();


    }

    public static class simParams extends ParameterBase{
        public int numAtoms = 500;
        public int numSteps = 200000;
    }


}

package etomica.modules.ensembles;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

public class LJMC extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono species;
    public final IBox box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveAtom mcMoveAtom;
    public final MCMoveVolume mcMoveVolume;
    public final MCMoveInsertDelete mcMoveID;

    public LJMC(Space _space) {
        super(_space);
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, 2.5, space);
        potentialMaster.setCellRange(2);
        int N = 182;  //number of atoms
        
        //controller and integrator
	    integrator = new IntegratorMC(potentialMaster, random, 1.0);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);//index 1
	    species.setIsDynamic(true);
        addSpecies(species);
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential,2.5);
	    potentialMaster.addPotential(p2Truncated, new IAtomType[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(space);
        addBox(box);
        IVectorMutable dim = space.makeVector();
        dim.E(space.D() == 2 ? 15 : 10);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);
        integrator.setBox(box);
        potentialMaster.getNbrCellManager(box).assignCellAll();
        
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
        
        mcMoveVolume = new MCMoveVolume(potentialMaster, random, space, 1) {
            public void acceptNotify() {
                if (vNew < Math.pow(6, space.D())) {
                    rejectNotify();
                    uNew = uOld;
                    return;
                }
                super.acceptNotify();
            }
        };
        
        mcMoveID = new MCMoveInsertDelete(potentialMaster, random, space) {
            public void acceptNotify() {
                if (moleculeList.getMoleculeCount()>999 && insert) {
                    rejectNotify();
                    uNew = 0;
                    return;
                }
                super.acceptNotify();
            }
        };
        mcMoveID.setSpecies(species);
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        System.out.println("initial energy "+meterPE.getDataAsScalar());
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        LJMC sim = new LJMC(space);
        sim.getController().actionPerformed();
    }
}

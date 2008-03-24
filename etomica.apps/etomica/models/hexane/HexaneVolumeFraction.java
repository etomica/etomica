package etomica.models.hexane;

/**
 * Class to calculate the volume of a single hexane molecule
 */
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;

public class HexaneVolumeFraction extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public IBox box;

    public BoundaryRectangularPeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
       
    public HexaneVolumeFraction(Space _space, int xCells, int yCells, int zCells) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(_space, false);
        PotentialMaster potentialMaster = new PotentialMaster(_space);
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = 6;
        primitive = new PrimitiveHexane(_space);
        lattice = new BravaisLattice(primitive);


        SpeciesHexane species = new SpeciesHexane(this, _space);
        getSpeciesManager().addSpecies(species);
        int[] nCells = new int[]{xCells, yCells, zCells};

//      int[] nCells = new int[]{10, 10, 10};
        bdry = new BoundaryRectangularPeriodic(getRandom(), _space);
        box = new Box(bdry, _space);
        addBox(box);
        box.setNMolecules(species, 1);
//        IVector tem = _space.makeVector();
//        tem.E(10.0);
//        System.out.println(tem);
//        box.setDimensions(tem);

         //Initialize the positions of the atoms.
        coordinateDefinition = new CoordinateDefinitionHexane(box, primitive, species, _space);
        coordinateDefinition.initializeCoordinates(nCells);
       
    }
    /**
     * @param args
     */
    public static void main(String[] args) {
        boolean graphic = true;
        int numberOfTests = 1000000;
        
        HexaneVolumeFraction sim = new HexaneVolumeFraction(Space3D.getInstance(), 2,2,2);
        
        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space);
            simGraphic.makeAndDisplayFrame();
        } else {
            int overlaps = 0;
            long time = System.currentTimeMillis();
            System.out.println(time);
            long time1;
            long time2;
            
            IVector rand = sim.space.makeVector();
            AtomIteratorLeafAtoms ail = new AtomIteratorLeafAtoms(sim.box);
            ail.reset();
            AtomLeaf atom = new AtomLeaf(sim.getSpace());
            
            time1 = System.currentTimeMillis();
            for(int count = 0; count < numberOfTests; count++){
                rand = ((BoundaryDeformablePeriodic)sim.box.getBoundary()).randomPosition();
                ail.reset();
                
                for(int atomNumber = 0; atomNumber < 512; atomNumber++){
                    atom = (AtomLeaf)ail.nextAtom();
                    
                }
            }
            time2 = System.currentTimeMillis();
            System.out.println("start  " + time);
            System.out.println("simulation  " + time1);
            System.out.println("data collection  " + time2);
        }
    }
}

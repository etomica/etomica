package etomica.models.hexane;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.lattice.BravaisLattice;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;

public class HexaneVolumeFraction extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public IBox box;

    public BoundaryDeformablePeriodic bdry;
    public BravaisLattice lattice;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
       
    public HexaneVolumeFraction(Space space, int xCells, int yCells, int zCells) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(space, false);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = xCells * yCells * zCells * chainLength;
        primitive = new PrimitiveHexane(space);
        // close packed density is 0.4165783882178116
        // Monson reports data for 0.373773507616 and 0.389566754417
        lattice = new BravaisLattice(primitive);


        SpeciesHexane species = new SpeciesHexane(this);
        getSpeciesManager().addSpecies(species);
        int[] nCells = new int[]{xCells, yCells, zCells};
        bdry = new BoundaryDeformableLattice(primitive, getRandom(), nCells);
        box = new Box(bdry, space);
        addBox(box);
        box.setNMolecules(species, xCells * yCells * zCells);
//        config.initializeCoordinates(box);

         //Initialize the positions of the atoms.
        coordinateDefinition = new CoordinateDefinitionHexane(box, primitive, species, space);
        coordinateDefinition.initializeCoordinates(nCells);
       
    }
    /**
     * @param args
     */
    public static void main(String[] args) {
        int xLng = 1;
        int yLng = 1;
        int zLng = 1;
        int numberOfTests = 1000000;
        
        HexaneVolumeFraction sim = new HexaneVolumeFraction(
                Space3D.getInstance(), xLng, yLng, zLng);
        
        
        int overlaps = 0;
        
        IVector rand = sim.space.makeVector();
        AtomIteratorLeafAtoms ail = new AtomIteratorLeafAtoms(sim.box);
        ail.reset();
        AtomLeaf atom = new AtomLeaf(sim.getSpace());
        
        
        for(int count = 0; count < numberOfTests; count++){
            rand = ((BoundaryDeformablePeriodic)sim.box.getBoundary()).randomPosition();
            ail.reset();
            
            for(int atomNumber = 0; atomNumber < 512; atomNumber++){
                atom = (AtomLeaf)ail.nextAtom();
                
            }
            
            
        }
            

    }

}

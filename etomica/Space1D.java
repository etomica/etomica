package etomica;

import etomica.atom.AtomTreeNodeGroup;
import etomica.space1d.*;
import etomica.statmech.MaxwellBoltzmann;

//CoordinateGroup is not updated to the same structure as used in Space2D and Space3D
//centralImage not updated to molecule form as in Space3D

 /* History of changes
  * 09/01/02 (DAK) added accelerateTo method to Coordinate
  *                changed CoordinateGroup.randomizeMomentum to not enforce zero COM momentum
  * 09/05/02 (DAK) fixed error in accelerateTo (still probably does not do what one expects
  *                if accelerating to nonzero momentum).
  * 07/10/03 (DAK) added resetV method to CoordinatePair
  * 08/27/03 (DAK) added isZero method to Vector
  * 08/29/03 (DAK) implemented centralImage(Space.Coordinate) in Boundary
  * 12/09/03 (DAK) changed setRandomSphere to give point on surface of sphere
  * (just at +/- 1 instead of between +/- 1); added setRandomInSphere method
  * 01/22/04 (DAK) added positionCOM and translateCOMTo to CoordinateGroup;
  * redefined position() in CoordinateGroup to be first-atom position (as it has
  * been for Space3D for some time now).
  */
public class Space1D extends Space implements EtomicaElement {
    
    public static final int D = 1;
    public static int drawingHeight = 10;  //height for drawing to 2D image
    public final int D() {return D;}
    public final int powerD(int n) {return n;}
    public final double powerD(double a) {return a;}
    public static final Vector ORIGIN = new Vector();
    public final etomica.space.Vector origin() {return ORIGIN;}
    public static final Space1D INSTANCE = new Space1D();
    
    public Space1D() {super(1);}
    
    public double sphereVolume(double r) {return 2.0*r;}  //volume of a sphere of radius r
    public double sphereArea(double r) {return 2.0;}      //surface area of sphere of radius r (used for differential shell volume)
    public etomica.space.Vector makeVector() {return new Vector();}
    public etomica.space.Orientation makeOrientation() {System.out.println("Orientation class not implemented in 1D"); return null;}
    public etomica.space.Tensor makeTensor() {return new Tensor();}
    public etomica.space.Tensor makeRotationTensor() {return new RotationTensor();}
    public etomica.space.Coordinate makeCoordinate(Atom a) {
        if(a.node instanceof AtomTreeNodeGroup) return new CoordinateGroup(a);
//        else if(a.type instanceof AtomType.Rotator) return new OrientedCoordinate(a);
        return new Coordinate(a);
    }
    public etomica.space.CoordinatePair makeCoordinatePair() {return new CoordinatePair();}

    public etomica.space.Boundary.Type[] boundaryTypes() {return Boundary.TYPES;}
    public etomica.space.Boundary makeBoundary() {return makeBoundary(Boundary.PERIODIC_SQUARE);}  //default
    public etomica.space.Boundary makeBoundary(etomica.space.Boundary.Type t) {
        if(t == Boundary.NONE) {return new BoundaryNone();}
        else if(t == Boundary.PERIODIC_SQUARE) {return new BoundaryPeriodicSquare();}
        else return null;
    }
    
    public int[] makeArrayD(int i) {return new int[] {i};}
    public double[] makeArrayD(double d) {return new double[] {d};}
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("One-dimensional space");
        return info;
    }

    public static final double r2(Vector u1, Vector u2, Boundary b) {
        Vector.WORK.x = u1.x - u2.x;
        b.nearestImage(Vector.WORK);
        return Vector.WORK.x*Vector.WORK.x;
    }
    
    
//    public static void main(String[] args) {
//        Default.ATOM_SIZE = 1.0;
//        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space1D());
//  //    setIteratorFactory(new IteratorFactoryCell(this));
//        AtomIteratorNeighbor nbrIterator = new AtomIteratorNeighbor();
//	    IntegratorMC integrator = new IntegratorMC(sim);
//	    MCMoveAtom mcmoveAtom = new MCMoveAtom(integrator);
//	    Species species = new SpeciesSpheresMono(sim);
//	    
//	    species.setNMolecules(3);
//	    final Phase phase = new Phase(sim);
//	    Potential2 potential = new P2HardSphere();
////	    Potential2 potentialOrder = new P2XOrder();
////	    potential.setSpecies(species,species);
////	    potentialOrder.setSpecies(species,species);
//	    Controller controller = new Controller(sim);
//	    etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(sim);
//        display.setColorScheme(new etomica.graphics.ColorSchemeRandom());
//		sim.elementCoordinator.go();
//		
//		AtomPairIterator api = new ApiGeneral(sim.space, new AtomIteratorList(), nbrIterator);
//		potential.setIterator(api);
//		
//		NeighborManager.Criterion criterion = new NeighborManager.Criterion() {
//		    public boolean areNeighbors(Atom a1, Atom a2) {
//		        return Math.abs(a1.node.index()-a2.node.index()) == 1 ||
//		           (a1==phase.firstAtom() && a2==phase.lastAtom());
//		    }};
//		nbrIterator.setupNeighbors(phase.speciesMaster.atomList, criterion);
//        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
//    }
            
}//end of Space1D

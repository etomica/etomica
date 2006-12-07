package etomica.junit.space;

import junit.framework.TestCase;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeGroup;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.IndexIteratorSequential;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;

public class BoundaryDeformablePeriodicTest extends TestCase {

    public BoundaryDeformablePeriodicTest() {
//        space = Space2D.getInstance();
        space = Space3D.getInstance();
        boundary = new BoundaryDeformablePeriodic(space, 1.0);
//        Tensor2D deformation = new Tensor2D(new double[] {1.0, -0.5, 
//                                                          0.8, 1.0});
//        Tensor2D deformation = new Tensor2D(new double[] {1.0, 12.0,
//      -.20, 1.0});
        Tensor3D deformation = new Tensor3D(new double[] {1.0, 0.10, 0.1,
                                                         -0.2, 1.0, .15,
                                                         -0.1, 0.5, 1.0});
        deformation.TE(30);
        boundary.deform(deformation);
        //System.out.println("in constructor, boundary="+boundary.getDimensions().toString());
        iMax = 50;
        dr = space.makeVector();
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        drStep = space.makeVector();
        imageIndexIterator = new IndexIteratorSequential(space.D());
        imageIndexIterator.setSize(3);
    }
    /*
     * Test method for 'etomica.space.BoundaryDeformablePeriodic.nearestImage(Vector)'
     */
    public void testNearestImage() {
        if(interactive) display.setPixelUnit(new Pixel(2));
        edgeVectors = space.makeVectorArray(space.D());  
        boundary.boundaryTensor().assignTo(edgeVectors);
        positionIndexIterator = new IndexIteratorSequential(space.D());
        positionIndexIterator.setSize(iMax);
        positionIndexIterator.reset();
        while(positionIndexIterator.hasNext()) {
            int[] index = positionIndexIterator.next();
            //index = new int[] {10, 8, 35};
            //System.out.println(Arrays.toString(index));
            dr.E(index);
            dr.TE(2.0/(double)iMax);
            dr.PE(-(1-1./iMax));
            dr1.E(dr);
            dr.transform(boundary.boundaryTensor());
            //System.out.println("dots: "+dr.dot(edgeVectors[0])/edgeVectors[0].squared()
            //                      +", "+dr.dot(edgeVectors[1])/edgeVectors[1].squared());
            dr1.E(dr);
            boundary.nearestImage(dr1);
            dr2.E(bruteForceNearestImage(dr));
            double delta = dr2.squared() - dr1.squared();
            if(delta < 1.e-10 && delta > -1.e-10) delta = 0.0;
            //           System.out.println(dr+" "+dr1+" "+dr2);
            if(interactive) {
                display.repaint();
                //System.out.println(dr+" "+dr1+" "+dr3);
                if(delta != 0.0) {
                    System.out.println(dr+" "+dr1+" "+dr2);
                    System.out.println("error"+dr.squared()+" "+dr1.squared()+" "+dr2.squared());
                    break;
                }
//                try {
//                    Thread.sleep(20);
//                } catch (InterruptedException e) {
//                    e.printStackTrace();
//                }
                //System.out.println();
            }
            assertTrue(delta == 0.0);
        }
    }
    
    public Vector bruteForceNearestImage(Vector dr) {
        double dr2Min = Double.MAX_VALUE;
        imageIndexIterator.reset();
        Vector drMin = (Vector)dr.clone();
        while(imageIndexIterator.hasNext()) {
            int[] idx = imageIndexIterator.next();
            drStep.E(dr);
            for(int i=0; i<idx.length; i++) {
                drStep.PEa1Tv1(idx[i]-1, edgeVectors[i]);
            }
            double dr2 = drStep.squared();
            if(dr2 < dr2Min) {
                dr2Min = dr2;
                drMin.E(drStep);
            }
        }
        return drMin;
    }
    
    public static SimulationGraphic makeDisplay(BoundaryDeformablePeriodicTest test) {
        Simulation sim = new Simulation(test.space);
        Phase phase = new Phase(sim);
        phase.setBoundary(test.boundary);
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        phase.getAgent(species).setNMolecules(3);
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DisplayPhase display = new DisplayPhase(phase, sim.getDefaults().pixelUnit);
        simGraphic.add(display);
        simGraphic.makeAndDisplayFrame();
        return simGraphic;
    }
    
    public static void main(String[] args) {
//        IndexIteratorSequential indexIterator = new IndexIteratorSequential(2);
//        int iMax = 10;
//        indexIterator.setSize(iMax);
//        indexIterator.reset();
//        Vector dr = Space2D.getInstance().makeVector();
//        while(indexIterator.hasNext()) {
//            int[] index = indexIterator.next();
//            dr.E(index);
//            dr.TE(2.0/(double)iMax);
//            dr.PE(-(1-1./iMax));
//           System.out.println(dr);
//        }
        BoundaryDeformablePeriodicTest test = new BoundaryDeformablePeriodicTest();
        test.simGraphic = makeDisplay(test);
        test.sim = test.simGraphic.getSimulation();
        test.atom0 = (AtomLeaf)((AtomTreeNodeGroup)test.sim.speciesRoot.getNode()).getDescendant(new int[] {0, 0, 0});
        test.atom1 = (AtomLeaf)((AtomTreeNodeGroup)test.sim.speciesRoot.getNode()).getDescendant(new int[] {0, 0, 1});
        test.atom2 = (AtomLeaf)((AtomTreeNodeGroup)test.sim.speciesRoot.getNode()).getDescendant(new int[] {0, 0, 2});
        test.display = ((DisplayPhase)test.simGraphic.displayList().getFirst());
        test.interactive = true;
        test.testNearestImage();
    }
    
    BoundaryDeformablePeriodic boundary;
    Space space;
    IndexIteratorSequential positionIndexIterator, imageIndexIterator;
    int iMax;
    Vector dr, dr1, dr2, drStep;
    Vector[] edgeVectors;
    Simulation sim;
    SimulationGraphic simGraphic;
    DisplayPhase display;
    AtomLeaf atom0, atom1, atom2;
    private boolean interactive = false;
}

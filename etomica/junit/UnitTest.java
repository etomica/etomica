package etomica.junit;

import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesRoot;
import etomica.SpeciesSpheres;
import etomica.SpeciesSpheresMono;
import etomica.SpeciesTree;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.space3d.Space3D;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Apr 28, 2005 by kofke
 */
public class UnitTest {

    public static boolean VERBOSE = false;
    
    /**
     * 
     */
    private UnitTest() {
        super();
        // TODO Auto-generated constructor stub
    }
    
    public static SpeciesRoot makeStandardSpeciesTree() {
        Space space = new Space3D();
        Simulation sim = new Simulation(space, new PotentialMaster(space), new int[] {1, 4, 4, 11, 6, 3, 3});
        SpeciesSpheres species0 = new SpeciesSpheres(sim, 3);
        species0.setNMolecules(5); //{5, 3}
        SpeciesSpheresMono species1 = new SpeciesSpheresMono(sim);
        species1.setNMolecules(10); //{10} 
        Species species2 = new SpeciesTree(sim, new int[] {5,4,3});
        species2.setNMolecules(3);//{3, 5, 4, 3}
        new Phase(sim);
        return sim.speciesRoot;
    }
    
    public static void main(String[] arg) {
        SpeciesRoot root = makeStandardSpeciesTree();
        AtomIteratorTree iterator = new AtomIteratorTree();
        iterator.setRoot(root);
        iterator.setDoAllNodes(true);
        iterator.reset();
        while(iterator.hasNext()) {
            System.out.println(iterator.next().toString());
        }
    }

}

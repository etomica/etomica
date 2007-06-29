package etomica.junit;

import junit.framework.TestCase;
import etomica.atom.AtomAddressManager;
import etomica.atom.AtomGroup;
import etomica.atom.AtomLeaf;
import etomica.atom.IAtom;
import etomica.atom.ISpeciesAgent;
import etomica.atom.SpeciesAgent;
import etomica.atom.AtomManager;
import etomica.box.Box;
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;


/**
 * Performs tests of AtomIndexManager, checking that methods to determine
 * if two atoms are in the same box, species, molecule, etc., function 
 * correctly.
 */
public class AtomIndexManagerTest extends TestCase {

    /*
     * @see TestCase#setUp()
     */
    protected void setUp() throws Exception {
        super.setUp();
        ISimulation sim = new Simulation(Space2D.getInstance());
        SpeciesSpheresMono species0 = new SpeciesSpheresMono(sim);
        SpeciesSpheres species1 = new SpeciesSpheres(sim, 5);
        sim.getSpeciesManager().addSpecies(species0);
        sim.getSpeciesManager().addSpecies(species1);
        Box box0 = new Box(sim);
        Box box1 = new Box(sim);
        sim.addBox(box0);
        sim.addBox(box1);
        box0.getAgent(species0).setNMolecules(20);
        box0.getAgent(species1).setNMolecules(10);
        box1.getAgent(species0).setNMolecules(20);
        box1.getAgent(species1).setNMolecules(10);
        atoms = new IAtom[24];
        int i = 0;
//        atoms[i++] = master0 = box0.getSpeciesMaster();//0
//        atoms[i++] = master1 = box1.getSpeciesMaster();//1
        atoms[i++] = agent00 = (SpeciesAgent)box0.getAgent(species0);//2
        atoms[i++] = agent01 = (SpeciesAgent)box0.getAgent(species1);//3
        atoms[i++] = agent10 = (SpeciesAgent)box1.getAgent(species0);//4
        atoms[i++] = agent11 = (SpeciesAgent)box1.getAgent(species1);//5
        AtomGroup[][] node = new AtomGroup[2][2];
        node[0][0] = agent00;
        node[0][1] = agent01;
        node[1][0] = agent10;
        node[1][1] = agent11;
        boxIndex =    new int[] {0,0,0,0,1,1,1,1, 0,0,0,0,0,0,1,1,1,1,1,1};
        speciesIndex =  new int[] {0,0,1,1,0,0,1,1, 1,1,1,1,1,1,1,1,1,1,1,1};
        moleculeIndex = new int[] {0,5,0,3,0,5,0,4, 0,0,3,3,4,4,0,0,3,3,4,4};
        atomIndex =     new int[] {0,5,0,3,0,5,0,4, 1,4,0,4,0,1,1,4,0,4,0,1};
        i0 = i;
        int del = 8;
        for(int j=0; j<del; j++) {
            atoms[i++] = node[boxIndex[j]][speciesIndex[j]].getChildList().getAtom(moleculeIndex[j]);//7
        }
        for(int j=del; j<boxIndex.length; j++) {
            atoms[i++] = node[boxIndex[j]][speciesIndex[j]].getDescendant(new int[] {moleculeIndex[j], atomIndex[j]});
        }
        atom = new AtomLeaf(Space3D.getInstance());
    }

    public void testGetOrdinal() {
        for(int i=0; i<atomIndex.length; i++) {
//          System.out.println("AtomIndex: "+atoms[i+i0].node.getOrdinal()+" "+(1+atomIndex[i]));
//          System.out.println("Bits: "+Integer.toBinaryString(atoms[i+i0].node.index()));
//          System.out.println(atoms[i+i0].toString());
          assertEquals(atoms[i+i0].getIndex(), atomIndex[i]);
      }
    }

    private static AtomManager getSpeciesMaster(IAtom atom) {
        IAtom speciesAgent = atom;
        while (!(speciesAgent instanceof ISpeciesAgent)) {
            speciesAgent = speciesAgent.getParentGroup();
        }
        return ((ISpeciesAgent)speciesAgent).getAtomManager();
    }
        
    
    public void testAncestry() {
        for(int i=0; i<atoms.length; i++) {
            assertFalse(atoms[i].isDescendedFrom(atom));
            assertFalse(atom.isDescendedFrom(atoms[i]));
            AtomManager iSpeciesMaster = getSpeciesMaster(atoms[i]);
            for(int j=0; j<atoms.length; j++) {
                if (iSpeciesMaster != getSpeciesMaster(atoms[j])) {
                    // different boxs, so skip
                    continue;
                }
//                System.out.println(i+" "+j);
                if(isDescendedFrom(atoms[i],atoms[j])) {
                    assertTrue(atoms[i].isDescendedFrom(atoms[j]));
                } else {
                    if (atoms[i].isDescendedFrom(atoms[j])) {
                        System.out.println("atoms i: "+Integer.toBinaryString(atoms[i].getAddress()));
                        System.out.println("atoms j: "+Integer.toBinaryString(atoms[j].getAddress()));
                        atoms[i].isDescendedFrom(atoms[j]);
                        System.out.println(atoms[i]+" "+atoms[j]);
                    }
                    assertFalse(atoms[i].isDescendedFrom(atoms[j]));
                }
            }
        }
    }
    private boolean isDescendedFrom(IAtom a1, IAtom a2) {
        if (a1 == null) return false;
        if(a1.getType().getDepth() < a2.getType().getDepth()) return false;
        else if(a1 == a2) return true;
        else return isDescendedFrom(a1.getParentGroup(), a2);
    }

    
    public void testIsDescendedFrom() {
        for(int i=0; i<atoms.length; i++) {
            assertFalse(atoms[i].getType().isDescendedFrom(atom.getType()));
            assertFalse(atom.getType().isDescendedFrom(atoms[i].getType()));
            AtomManager iSpeciesMaster = getSpeciesMaster(atoms[i]);
            for(int j=0; j<atoms.length; j++) {
                if (iSpeciesMaster != getSpeciesMaster(atoms[j])) {
                    // different boxs, so skip
                    continue;
                }
                boolean is = atoms[i].getType().isDescendedFrom(atoms[j].getType());
                if(typeIsDescendedFrom(atoms[i],atoms[j])) {
                    assertTrue(is);
                } else {
                    if(is) {
                        System.out.println("isDescendedFrom "+i+" "+j+" "+atoms[i]+" "+atoms[j]+" "+Integer.toBinaryString(atoms[i].getType().getAddress())+" "+Integer.toBinaryString(atoms[j].getType().getAddress()));
                        typeIsDescendedFrom(atoms[i],atoms[j]);
                    }
                    assertFalse(is);
                }
            }
        }
    }
    private boolean typeIsDescendedFrom(IAtom a1, IAtom a2) {
        if (a1 == null) return false;
        if(a1.getType().getDepth() < a2.getType().getDepth()) return false;
        else if(a1.getType() == a2.getType()) return true;
        else return typeIsDescendedFrom(a1.getParentGroup(), a2);
    }
    
    public void testSameMolecule() {
        int trueCount = 0;
        int falseCount = 0;
        int undefinedCount = 0;
        for(int i=i0; i<atoms.length; i++) {
            if (atoms[i].getType().getDepth() < AtomAddressManager.MOLECULE_DEPTH) {
                continue;
            }
            if(atoms[i].inSameMolecule(atom)) System.out.println(i+" "+atoms[i]+" "+atom);
            assertFalse(atoms[i].inSameMolecule(atom));
            assertFalse(atom.inSameMolecule(atoms[i]));
            IAtom moleculeA = atoms[i];
            while (!(moleculeA.getParentGroup() instanceof ISpeciesAgent)) {
                moleculeA = moleculeA.getParentGroup();
            }
            AtomManager iSpeciesMaster = getSpeciesMaster(atoms[i]);
            for(int j=i0; j<atoms.length; j++) {
                if (atoms[j].getType().getDepth() < AtomAddressManager.MOLECULE_DEPTH) {
                    continue;
                }
                if (iSpeciesMaster != getSpeciesMaster(atoms[j])) {
                    // different boxs, so skip
                    continue;
                }
                IAtom moleculeB = atoms[j];
                while (!(moleculeB.getParentGroup() instanceof ISpeciesAgent)) {
                    moleculeB = moleculeB.getParentGroup();
                }
                boolean inSameMolecule = atoms[i].inSameMolecule(atoms[j]);
                if(moleculeA == null || moleculeB == null) {
                    if(inSameMolecule) System.out.println(inSameMolecule+" "+i+" "+j+" "+atoms[i]+" "+atoms[j]);
                    assertFalse(inSameMolecule);
                    undefinedCount++;
                } else if(moleculeA == moleculeB) {
                    if(!inSameMolecule) System.out.println(inSameMolecule+" "+i+" "+j+" "+atoms[i]+" "+atoms[j]);
                    assertTrue(inSameMolecule);
                    trueCount++;
                } else {
                    if(inSameMolecule) System.out.println(inSameMolecule+" "+i+" "+j+" "+atoms[i]+" "+atoms[j]);
                    assertFalse(inSameMolecule);
                    falseCount++;
                }
            }
        }
        if(UnitTestUtil.VERBOSE) System.out.println("AtomIndexManagerTest(Molecule): trueCount = "+trueCount+"; falseCount = "+falseCount+"; undefinedCount = "+undefinedCount);
    }

    IAtom atom;
    SpeciesAgent agent00, agent01, agent10, agent11;
    AtomManager master0, master1;
    IAtom[] atoms;
    int[] moleculeIndex, boxIndex, speciesIndex, atomIndex;
    int i0;
}

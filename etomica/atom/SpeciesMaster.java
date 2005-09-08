package etomica.atom;

import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

/**
 * Coordinator of all species agents in a phase. Parent is SpeciesRoot, and
 * children are SpeciesAgent instances. All instances of SpeciesMaster in a
 * given simulation share the same AtomType.
 * 
 * @author David Kofke and Andrew Schultz
 */

public final class SpeciesMaster extends Atom {

    protected int moleculeCount;
    public final static int SPECIES_TAB = AtomLinker.Tab.requestTabType();
    //reference to phase is kept in node

    /**
     * Tabbed list of leaf atoms in phase, suitable for iteration via an
     * AtomIteratorTabbedList.
     */
    public final AtomList atomList = new AtomListTabbed();

    public SpeciesMaster(Simulation sim, Phase p) {
        super(null, sim.speciesRoot.getChildType(), new NodeFactory(p),
                AtomLinker.FACTORY);
    }

    public SpeciesAgent firstSpecies() {
        return (SpeciesAgent) ((AtomTreeNodeGroup)node).childList.getFirst();
    }

    public SpeciesAgent lastSpecies() {
        return (SpeciesAgent) ((AtomTreeNodeGroup)node).childList.getLast();
    }

    //    public int atomCount() {return atomList.size();}//or could use
    // node.leafAtomCount()
    public int moleculeCount() {
        return moleculeCount;
    }

    public String signature() {
        Phase phase = node.parentPhase();
        return (phase != null) ? phase.getName()
                : "SpeciesMaster without phase";
    }

    private static final class MasterAtomTreeNode extends AtomTreeNodeGroup {

        MasterAtomTreeNode(Phase parentPhase, Atom atom) {
            super(atom);
            speciesMaster = (SpeciesMaster) atom;
            this.parentPhase = parentPhase;
            leafIterator.setAsLeafIterator();
        }

        public Phase parentPhase() {
            return parentPhase;
        }

        public Species parentSpecies() {
            throw new RuntimeException(
                    "Error:  Unexpected call to parentSpecies in SpeciesMaster");
        }

        public SpeciesAgent parentSpeciesAgent() {
            throw new RuntimeException(
                    "Error:  Unexpected call to parentSpeciesAgent in SpeciesMaster");
        }

        /**
         * Throws a RuntimeException, because a species master is not contained
         * within a molecule.
         */
        public final Atom parentMolecule() {
            throw new RuntimeException(
                    "Error:  Unexpected call to parentMolecule in SpeciesMaster");
        }

        /**
         * Returns true, because children are SpeciesAgent instances.
         */
        public final boolean childrenAreGroups() {
            return true;
        }

        //does not pass notification up to parent (root)
        public void addAtomNotify(Atom newAtom) {
            if (newAtom.node.parentGroup() instanceof SpeciesAgent) {
                speciesMaster.moleculeCount++;
            } else if (newAtom instanceof SpeciesAgent) {
                speciesMaster.moleculeCount += ((SpeciesAgent) newAtom)
                        .moleculeCount();
                AtomLinker.Tab newTab = AtomLinker.newTab(
                        speciesMaster.atomList, SPECIES_TAB);
                speciesMaster.atomList.add(newTab);
                ((SpeciesAgent) newAtom).firstLeafAtomTab = newTab;
            }
            AtomLinker.Tab nextTab = newAtom.node.parentSpeciesAgent().firstLeafAtomTab.nextTab;

            leafAtomCount += newAtom.node.leafAtomCount();
            if (newAtom.node.isLeaf()) {
                speciesMaster.atomList.addBefore(
                        ((AtomTreeNodeLeaf) newAtom.node).leafLinker, nextTab);
            } else {
                leafIterator.setRoot(newAtom);
                leafIterator.reset();
                while (leafIterator.hasNext()) {
                    speciesMaster.atomList
                            .addBefore(((AtomTreeNodeLeaf) leafIterator
                                    .nextAtom().node).leafLinker, nextTab);
                }
            }
       }

        //updating of leaf atomList may not be efficient enough for repeated
        // use, but is probably ok
        public void removeAtomNotify(Atom oldAtom) {
            if (oldAtom.node.parentGroup() instanceof SpeciesAgent) {
                speciesMaster.moleculeCount--;
            } else if (oldAtom instanceof SpeciesAgent) {
                speciesMaster.moleculeCount -= ((SpeciesAgent) oldAtom)
                        .moleculeCount();
            }
            if (oldAtom.node.isLeaf()) {
                leafAtomCount--;
                speciesMaster.atomList.remove(((AtomTreeNodeLeaf) oldAtom.node).leafLinker);

            } else {
                leafAtomCount -= oldAtom.node.leafAtomCount();
                leafIterator.setRoot(oldAtom);
                leafIterator.reset();
                while (leafIterator.hasNext()) {
                    speciesMaster.atomList
                            .remove(((AtomTreeNodeLeaf) leafIterator.nextAtom().node).leafLinker);
                }
            }
        }

        private final Phase parentPhase;
        private final SpeciesMaster speciesMaster;
        private final AtomIteratorTree leafIterator = new AtomIteratorTree();
    } //end of MasterAtomTreeNode

    private static final class NodeFactory implements AtomTreeNodeFactory, java.io.Serializable {

        Phase phase;

        NodeFactory(Phase p) {
            phase = p;
        }

        public AtomTreeNode makeNode(Atom atom) {
            return new MasterAtomTreeNode(phase, atom);
        }
    }

    /**
     * non-graphic main method to test handling of leaf atom list.
     */
    public static void main(String args[]) {

        Simulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono(sim);
        Species species1 = new SpeciesSpheres(sim, 3);
        Species species0 = new SpeciesSpheres(sim, 2);
        species0.setNMolecules(4);
        species1.setNMolecules(2);
        species2.setNMolecules(2);
        Phase phase = new Phase(sim);
        //        sim.elementCoordinator.go();

        AtomIteratorLeafAtoms leafIterator = new AtomIteratorLeafAtoms();
        leafIterator.setPhase(phase);
        leafIterator.reset();
        while (leafIterator.hasNext())
            System.out.println(leafIterator.next().toString());
        System.out.println();
    }//end of main

}//end of SpeciesMaster

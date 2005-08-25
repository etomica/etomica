package etomica.species;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryTree;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactory;
import etomica.space.CoordinateFactorySphere;

/**
 * Species in which molecules are formed as an arbitrarily shaped tree.
 * 
 * @author David Kofke
 */

public class SpeciesTree extends Species implements EtomicaElement {

    /**
     * Constructs with nA = {1}, such that each molecule is a group
     * containing just one atom (which is not the same as SpeciesSpheresMono,
     * for which each molecule is a single atom, not organized under a group).
     */
    public SpeciesTree(Simulation sim) {
        this(sim, new int[] {1});
    }

    /**
     * Constructor specifing tree structure through array of integers.
     * Each element of array indicates the number of atoms at the corresponding
     * level.  For example, nA = {2,4} will define a species in which each
     * molecule has 2 subgroups, each with 4 atoms (such as ethane, which
     * can be organized as CH3 + CH3)
     */
    public SpeciesTree(Simulation sim, int[] nA) {
        this(sim, sim.potentialMaster.sequencerFactory(), nA);
    }
    
    public SpeciesTree(Simulation sim, AtomSequencerFactory seqFactory, int[] nA) {
        this(sim, seqFactory, nA, Species.makeAgentType(sim));
    }
    
    //TODO extend to permit specification of Conformation[], perhaps AtomSequencerFactory[]
    private SpeciesTree(Simulation sim, AtomSequencerFactory seqFactory, 
            int[] nA, AtomTypeGroup agentType) {
        super(sim, new AtomFactoryTree(sim.space, seqFactory, agentType, nA), agentType);
        //getLeafType will return the an AtomTypeGroup because leaf factory is not yet set
        AtomTypeSphere atomType = new AtomTypeSphere((AtomTypeGroup)((AtomFactoryTree)factory).getLeafType(), sim.getDefaults().atomMass, sim.getDefaults().atomSize);
        CoordinateFactory coordFactory = new CoordinateFactorySphere(sim);
        ((AtomFactoryTree)factory).setLeafFactory(new AtomFactoryMono(coordFactory, atomType, seqFactory));
        factory.setSpecies(this);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules formed as an arbitrarily specified tree");
        return info;
    }
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{Simulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{new Integer(1)});
    }
}



/*
 * Created on May 18, 2005
 */
package etomica.models.hexane;

import etomica.Atom;
import etomica.AtomIndexManager;
import etomica.AtomType;
//import etomica.Conformation;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.action.AtomActionAdapter;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeSphere;
import etomica.units.Dimension;

/**
 * Species used to create hexane molecules per Dr. Monson's data. Hydrogen
 * molecules are ignored.
 * 
 * @author nancycribbin
 */

public class SpeciesHexane extends Species implements EtomicaElement {

    private double mass;

    private final AtomTypeSphere atomType;

    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryHomo makeFactory(Space space,
            AtomSequencerFactory seqFactory, AtomIndexManager indexManager,
            int na, ConformationHexane conform) {
        AtomType type = new AtomTypeSphere(indexManager.makeChildManager(),
                Default.ATOM_MASS, Default.ATOM_SIZE);
        AtomFactoryMono f = new AtomFactoryMono(space, type, seqFactory);
        return new AtomFactoryHomo(space, seqFactory, indexManager, f, na,
                conform);
    }

    public SpeciesHexane(Simulation sim) {
        this(sim, sim.potentialMaster.sequencerFactory(), 1);
    }

    public SpeciesHexane(Simulation sim, AtomSequencerFactory seqFactory, int nA) {
        this(sim, seqFactory, nA, new ConformationHexane(sim.space));
    }

    public SpeciesHexane(Simulation sim, AtomSequencerFactory seqFactory,
            int nA, ConformationHexane conformation) {
        this(sim, seqFactory, nA, conformation, Species.makeAgentType(sim));
    }

    private SpeciesHexane(Simulation sim, AtomSequencerFactory seqFactory,
            int nA, ConformationHexane conformation, AtomType agentType) {
        super(sim, makeFactory(sim.space, seqFactory, agentType
                .getIndexManager().makeChildManager(), nA, conformation),
                agentType);
        factory.setSpecies(this);
        atomType = (AtomTypeSphere) ((AtomFactoryHomo) factory).childFactory()
                .getType();
        nMolecules = Default.MOLECULE_COUNT;
        mass = atomType.getMass();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hexane chain species.");
        return info;
    }

    /**
     * The mass of each of the spheres that form a molecule.
     */
    public final double getMass() {
        return mass;
    }

    /**
     * Sets the mass of all spheres in each molecule to the given value.
     */
    public final void setMass(double m) {
        mass = m;
        atomType.setMass(m);
    }

    /**
     * @return Dimension.MASS
     */
    public Dimension getMassDimension() {
        return Dimension.MASS;
    }

    /**
     * The diameter of each of the spheres that form a molecule.
     */
    public final double getDiameter() {
        return atomType.diameter(null);
    }

    /**
     * Sets the diameter of all spheres in each molecule to the given value.
     */
    public void setDiameter(double d) {
        atomType.setDiameter(d);
    }

    /**
     * @return Dimension.LENGTH
     */
    public Dimension getDiameterDimension() {
        return Dimension.LENGTH;
    }

    /**
     * Sets the number of spheres in each molecule. Causes reconstruction of all
     * molecules of this species in all phases. No action is performed if the
     * given value equals the current value.
     * 
     * @param n
     *            new number of atoms in each molecule.
     */
    public void setAtomsPerMolecule(final int n) {
        if (n == getAtomsPerMolecule())
            return;
        ((AtomFactoryHomo) factory).setAtomsPerGroup(n);
        allAgents(new AtomActionAdapter() {
            public void actionPerformed(Atom a) {
                SpeciesAgent agent = (SpeciesAgent) a;
                agent.setNMolecules(agent.getNMolecules());
            }
        });
    }

    /**
     * @return the number of spheres in each molecule made by this species.
     */
    public int getAtomsPerMolecule() {
        return ((AtomFactoryHomo) factory).getAtomsPerGroup();
    }

}
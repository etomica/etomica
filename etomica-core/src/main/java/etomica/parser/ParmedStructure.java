package etomica.parser;


import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IVector;
import etomica.atom.AtomTypeLeaf;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresCustom;

import java.util.*;

/**
 * Class that contains the data from a <a href="https://github.com/ParmEd/ParmEd">ParmEd</a> {@code Structure},
 * and has methods to load this data into various <i>etomica</i> components.
 *
 * <p>
 * This class should not be instantiated directly, but instead created by the {@link ParmedParser}
 * static methods which invoke the Python library.
 * </p>
 *
 * @see ParmedParser
 */
public class ParmedStructure {
    private static final ObjectMapper mapper = new ObjectMapper();
    private final JsonNode root;

    // We're only ever using 3D space for now
    private static final Space SPACE = Space3D.getInstance();

    private Map<String, AtomTypeLeaf> atomTypes;
    private SpeciesSpheresCustom species;

    /**
     * rmin = 2^(1/6) * sigma
      */
    private static final double SIGMA_FACTOR = 1 / Math.pow(2, 1/6);

    ParmedStructure(JsonNode root) {
        this.root = root;
    }

    /**
     * Gets the (single, rectangular) box found in this Structure.
     *
     * @return a Box object with the same boundary coordinates as the box in the Structure.
     */
    public Box getBox() {
        JsonNode boxNode = root.get("_box");
        double[] boxGeometry;
        try {
            boxGeometry = mapper.treeToValue(boxNode.get(0), double[].class);
        } catch (JsonProcessingException e) {
            throw new RuntimeException("Error during JSON processing");
        }
        double[] boxCoordinates = Arrays.copyOf(boxGeometry, 3);

        Box box = new Box(SPACE);
        Boundary bound = new BoundaryRectangularPeriodic(SPACE, boxCoordinates);
        box.setBoundary(bound);
        return box;
    }

    public SpeciesSpheresCustom getSpecies() {
        if(!Objects.isNull(species)) {
            return this.species;
        }

        ensureAtomTypes();

        SpeciesSpheresCustom theSpecies = new SpeciesSpheresCustom(
                SPACE,
                atomTypes.values().toArray(new AtomTypeLeaf[]{ })
        );

        // LinkedHashMap keySet is guaranteed to be in insertion order
        List<String> atomTypeIndices = new ArrayList<>(atomTypes.keySet());
        List<Integer> speciesAtomTypes = new ArrayList<>();
        List<IVector> atomPositions = new ArrayList<>();
        JsonNode speciesAtoms = root.get("residues").get(0).get("atoms");

        for(JsonNode atomNode : speciesAtoms) {
            String atomType = atomNode.get("type").asText();
            speciesAtomTypes.add(atomTypeIndices.indexOf(atomType));

            atomPositions.add(new Vector3D(
                    atomNode.get("xx").asDouble(),
                    atomNode.get("xy").asDouble(),
                    atomNode.get("xz").asDouble()
            ));

        }

        // need to convert the list of Integers into array of ints
        theSpecies.setAtomTypes(speciesAtomTypes.stream().mapToInt(i -> i).toArray());

        theSpecies.setConformation(atomList -> {
            for(int i = 0; i < atomList.getAtomCount(); i++) {
                IVector atomVec = atomPositions.get(i);
                atomList.getAtom(i).getPosition().E(atomVec);
            }
        });

        this.species = theSpecies;
        return theSpecies;
    }

    /**
     * Gets the intermolecular {@link PotentialGroup} for the atom types found in this Structure. Constructs
     * all unique unordered pairs of atom types and creates {@link P2LennardJones} potentials using
     * <a href="https://en.wikipedia.org/wiki/Combining_rules">Lorentz-Berthelot combining rules</a>.
     *
     * @return the PotentialGroup for potentials between molecules in this Structure.
     */
    public PotentialGroup getIntermolecularPotential() {
        ensureAtomTypes();
        PotentialGroup potentialGroup = new PotentialGroup(2, SPACE);

        JsonNode atomTypesNode = root.get("parameterset").get("atom_types");
        List<JsonNode> atomTypesList = new ArrayList<>();
        atomTypesNode.elements().forEachRemaining(atomTypesList::add);

        for (int i = 0; i < atomTypesList.size(); i++) {
            for (int j = i; j < atomTypesList.size(); j++) {

                JsonNode typeNode1 = atomTypesList.get(i);
                JsonNode typeNode2 = atomTypesList.get(j);

                AtomTypeLeaf atomType1 = atomTypes.get(typeNode1.get("name").asText());
                AtomTypeLeaf atomType2 = atomTypes.get(typeNode2.get("name").asText());
                AtomTypeLeaf[] typePair = new AtomTypeLeaf[] { atomType1, atomType2 };

                double epsilon1 = typeNode1.get("epsilon").asDouble();
                double epsilon2 = typeNode2.get("epsilon").asDouble();
                // geometric mean
                double combinedEpsilon = Math.sqrt(epsilon1 * epsilon2);

                double sigma1 = typeNode1.get("rmin").asDouble() * SIGMA_FACTOR;
                double sigma2 = typeNode2.get("rmin").asDouble() * SIGMA_FACTOR;
                double combinedSigma = (sigma1 + sigma2) / 2;

                P2LennardJones potential = new P2LennardJones(SPACE, combinedSigma, combinedEpsilon);

                potentialGroup.addPotential(potential, typePair);
            }
        }

        return potentialGroup;
    }

    /**
     * Gets the list of all molecules in the Structure, initialized with their coordinates in space.
     *
     * @return a List containing all the molecules of a species in the structure.
     */
    public List<IMolecule> getMolecules() {
        List<IMolecule> moleculeList = new ArrayList<>();
        SpeciesSpheresCustom species = getSpecies();

        //TODO: extract by species, i.e. don't assume all molecules are same species
        for(JsonNode moleculeNode : root.get("residues")) {
            IMolecule molecule = species.makeMolecule();

            JsonNode atomsListNode = moleculeNode.get("atoms");
            IAtomList atomList = molecule.getChildList();
            for(int i = 0; i < atomList.getAtomCount(); i++) {
                JsonNode atomNode = atomsListNode.get(i);
                atomList.getAtom(i).getPosition().E(new Vector3D(
                        atomNode.get("xx").asDouble(),
                        atomNode.get("xy").asDouble(),
                        atomNode.get("xz").asDouble()
                ));

            }
            moleculeList.add(molecule);
        }
        return moleculeList;
    }

    /**
     * "Lazily" create the atomTypes field when it is first needed. If the field is null it creates the
     * map from the data, otherwise it does nothing.
     */
    private void ensureAtomTypes() {
        if(Objects.isNull(atomTypes)) {
            JsonNode atomTypesList = root.get("parameterset").get("atom_types");
            atomTypes = new LinkedHashMap<>();

            for(JsonNode atomTypeNode : atomTypesList) {
                String atomName = atomTypeNode.get("name").asText();
                double atomMass = atomTypeNode.get("mass").asDouble();

                ElementSimple element = new ElementSimple(atomName, atomMass);
                AtomTypeLeaf atomType = new AtomTypeLeaf(element);
                atomTypes.put(atomName, atomType);
            }
        }
        // else do nothing
    }
}

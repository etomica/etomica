package etomica.parser;


import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.api.IVector;
import etomica.atom.AtomTypeLeaf;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresCustom;

import java.util.*;

public class ParmedStructure {
    private static final ObjectMapper mapper = new ObjectMapper();
    private final JsonNode root;

    // We're only ever using 3D space for now
    private static final Space simSpace = Space3D.getInstance();

    ParmedStructure(JsonNode root) {
        this.root = root;
    }

    public Box getBox() throws JsonProcessingException {
        JsonNode boxNode = root.get("_box");
        double[] boxGeometry = mapper.treeToValue(boxNode.get(0), double[].class);
        double[] boxCoordinates = Arrays.copyOf(boxGeometry, 3);

        Box box = new Box(simSpace);
        Boundary bound = new BoundaryRectangularPeriodic(simSpace, boxCoordinates);
        box.setBoundary(bound);
        return box;
    }

    public SpeciesSpheresCustom getSpecies() {
        JsonNode atomTypesList = root.get("parameterset").get("atom_types");
        Map<String, AtomTypeLeaf> atomTypes = new LinkedHashMap<>();

        for(JsonNode atomTypeNode : atomTypesList) {
            String atomName = atomTypeNode.get("name").asText();
            double atomMass = atomTypeNode.get("mass").asDouble();

            ElementSimple element = new ElementSimple(atomName, atomMass);
            AtomTypeLeaf atomType = new AtomTypeLeaf(element);
            atomTypes.put(atomName, atomType);
        }

        SpeciesSpheresCustom theSpecies = new SpeciesSpheresCustom(
                simSpace,
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

        return theSpecies;
    }
}

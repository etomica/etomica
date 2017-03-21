package etomica.parser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import etomica.atom.AtomTypeLeaf;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresCustom;

/**
 * Created by alex on 3/16/17.
 */
public class ParmedParser {
    private static final ObjectMapper mapper = new ObjectMapper();

    public static void main(String[] args) {
        try {
            File topFile = getResourceFile("test.top");
            File groFile = getResourceFile("test.gro");


            ProcessBuilder pb = new ProcessBuilder(
                    "venv/bin/parmed_json",
                    topFile.getCanonicalPath(),
                    groFile.getCanonicalPath()
            );


            Process p = pb.start();

            JsonNode tree = mapper.readTree(p.getInputStream());

            p.waitFor();

//            System.out.println(tree);
            JsonNode boxNode = tree.get("_box");
            double[] b = mapper.treeToValue(boxNode.get(0), double[].class);
            Box box = new Box(Space3D.getInstance());
            Boundary bound = new BoundaryRectangularPeriodic(Space3D.getInstance(), Arrays.copyOf(b, 3));
            box.setBoundary(bound);
            
      
            JsonNode atomTypesList = tree.get("parameterset").get("atom_types");
            System.out.println(atomTypesList);
            Map<String, AtomTypeLeaf> atomTypes = new LinkedHashMap<>();
            List<String> atomTypeIndices = new ArrayList<>();
            
            for(JsonNode atomType : atomTypesList) {
            	System.out.println(atomType);
            	ElementSimple el = new ElementSimple(atomType.get("name").asText(), atomType.get("mass").asDouble());
            	AtomTypeLeaf at = new AtomTypeLeaf(el);
            	atomTypes.put(atomType.get("name").asText(), at);
            	atomTypeIndices.add(atomType.get("name").asText());
            }
            
            JsonNode speciesAtoms = tree.get("residues").get(0).get("atoms");
            SpeciesSpheresCustom theSpecies = new SpeciesSpheresCustom(Space3D.getInstance(), atomTypes.values().toArray(new AtomTypeLeaf[]{}));
            List<Integer> speciesAtomTypes = new ArrayList<>();
            for(JsonNode atom : speciesAtoms) {
            	speciesAtomTypes.add(atomTypeIndices.indexOf(atom.get("type").asText()));
            	
            }
            
            theSpecies.setAtomTypes(speciesAtomTypes.stream().mapToInt(i -> i).toArray());
            
          
            

        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static File getResourceFile(String filename) {
        ClassLoader classLoader = ParmedParser.class.getClassLoader();
        return new File(classLoader.getResource(filename).getFile());
    }
}

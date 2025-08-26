package etomica.atom;

import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import java.util.ArrayList;
import java.util.List;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;


public class buildMethane {

    public static ISpecies buildMethane(String confName){
        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        String fileName = confName + ".pdb";
        FileReader fileReader;
        Map<String, AtomType> typeMap = new HashMap<>();
        List<List<Integer>> connectivity = new ArrayList<>();
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open " + fileName + ", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line;

            while ((line = bufReader.readLine()) != null) {
                parseLine(line, speciesBuilder, typeMap);

                if (line.startsWith("CONECT")) {
                    String[] parts = line.trim().split("\\s+");
                    int atomNumber = Integer.parseInt(parts[1]);
                    List<Integer> currentAtomList = null;

                    for (List<Integer> atomList : connectivity) {
                        if (atomList.get(0) == atomNumber) {
                            currentAtomList = atomList;
                            break;
                        }
                    }

                    if (currentAtomList == null) {
                        currentAtomList = new ArrayList<>();
                        connectivity.add(currentAtomList);
                        currentAtomList.add(atomNumber);
                    }

                    for (int i = 2; i < parts.length; i++) {
                        currentAtomList.add(Integer.parseInt(parts[i]));
                    }
                }
            }
            fileReader.close();
        } catch (IOException e) {
            throw new RuntimeException("Problem reading from " + fileName + ", caught IOException: " + e.getMessage());
        }

        // Assuming the SpeciesBuilder is configured properly, create the Species object
        return speciesBuilder.build(); // Adjusted this line
    }

    private static void parseLine(String line, SpeciesBuilder speciesBuilder, Map<String, AtomType> typeMap) {
        // Implement this method to parse lines from the PDB file and configure the SpeciesBuilder accordingly
    }
}

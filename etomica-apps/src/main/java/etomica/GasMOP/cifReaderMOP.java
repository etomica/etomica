package etomica.GasMOP;

import etomica.atom.AtomType;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class cifReaderMOP {
/*
    public void readCIFFile(String confName, double multiplier){
        String fileName = confName+".cif";
        FileReader fileReader;
        double[] boxSize = new double[3];
        double[] boxAngle = new double[3];
        operations = parseLineReader(fileName, boxSize, boxAngle, positionCoordinatesMap, atomTypeMap);
        System.out.println("\n");
        System.out.println(boxSize);
        // System.out.println(Arrays.deepToString(operations.toArray()));
        System.out.println(atomTypeMap);
        // System.out.println(positionCoordinatesMap);
        ArrayList<double[]> newSet = new ArrayList<>();
        int i=0;
        for(String key : positionCoordinatesMap.keySet()){
            double[] vectArrayPosn = positionCoordinatesMap.get(key).toArray();
            newSet = getCoordinates(vectArrayPosn, operations, boxSize, multiplier);
            atomTypePositions.put(key, newSet);
            i++;
        }
      /*  for (Map.Entry<String, ArrayList<double[]>> entry : atomTypePositions.entrySet()) {
            String key = entry.getKey();
            ArrayList<double[]> positions = entry.getValue();
         //   System.out.println("Atom Type: " + key);
            /*for (double[] position : positions) {
                System.out.println(" Position: " + Arrays.toString(position));
            }
        }
    }


    public ISpecies speciesCIF (String confName, boolean isDynamic){
        AtomType typeNew;
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        readCIFFile(confName, 5);
        int i=0;
        for (String key: positionCoordinatesMap.keySet()){
            String atomType = atomTypeMap.get(key);
            Pattern pattern = Pattern.compile("([a-zA-Z]+)(\\d+)");
            Matcher matcher = pattern.matcher(atomType);
            if (matcher.find()) {
                String letters = matcher.group(1);
                String numbers = matcher.group(2);
                i= Integer.parseInt(numbers);
            }else {
                //  System.out.println("Pattern Not found");
            }
            AtomType newAtom = returnElement(key, isDynamic);
            System.out.println(atomType);
            if (typeMapNew.containsKey(atomType)) {
                typeNew = typeMapNew.get(atomType);
            } else {
                typeNew = newAtom;
                typeMapNew.put(atomType, typeNew);
            }
            ArrayList<double[]> atomPositions = atomTypePositions.get(key);
            //  System.out.println(atomPositions);
            for(int j=0; j<atomPositions.size(); j++){
                Vector position = Vector.of(atomPositions.get(j));
                //System.out.println(position);
                speciesBuilderNew.addAtom(typeNew, position,  "");
            }
        }
        System.out.println("Done Building");
        speciesCIF = speciesBuilderNew.setDynamic(isDynamic).build();
        return speciesCIF;
    }*/
}

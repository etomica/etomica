package etomica.potential.UFF;

import etomica.atom.AtomType;
import etomica.chem.elements.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class cifReader {
    static Map<String, ArrayList<double[]>> positionList = new HashMap<>();
    Map<String, Vector> positionCoordinatesMap = new HashMap<>();
    Map<String, String> atomTypeMap = new HashMap<>();
    Map<String, AtomType> elementReceiverMap = new HashMap<>();
    ISpecies speciesCIF;
    ArrayList<String[]> operations = new ArrayList<>();
    Map<String, ArrayList<double[]>> atomTypePositions = new HashMap<>();
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
    }
    Map<String, AtomType> typeMapNew = new HashMap<>();
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
        }*/
    }
    public AtomType returnElement(String element, boolean isInfinite){
        if(!isInfinite){
            if(element.equals("O1")){
                elementReceiverMap.put("O1", new AtomType(Oxygen.INSTANCE, "O1"));
            }else if (element.equals("O2")){
                elementReceiverMap.put("O2", new AtomType(Oxygen.INSTANCE, "O2"));
            } else if (element.equals("O3")) {
                elementReceiverMap.put("O3", new AtomType(Oxygen.INSTANCE, "O3"));
            } else if (element.equals("C1")) {
                elementReceiverMap.put("C1", new AtomType(Carbon.INSTANCE, "C1"));
            } else if (element.equals("C2")) {
                elementReceiverMap.put("C2", new AtomType(Carbon.INSTANCE, "C2"));
            }else if (element.equals("C3")) {
                elementReceiverMap.put("C3", new AtomType(Carbon.INSTANCE, "C3"));
            }else if (element.equals("C4")) {
                elementReceiverMap.put("C4", new AtomType(Carbon.INSTANCE, "C4"));
            }else if (element.equals("C5")) {
                elementReceiverMap.put("C5", new AtomType(Carbon.INSTANCE, "C5"));
            }else if (element.equals("C6")) {
                elementReceiverMap.put("C6", new AtomType(Carbon.INSTANCE, "C6"));
            }else if (element.equals("C7")) {
                elementReceiverMap.put("C7", new AtomType(Carbon.INSTANCE, "C7"));
            }else if (element.equals("C8")) {
                elementReceiverMap.put("C8", new AtomType(Carbon.INSTANCE, "C8"));
            }else if (element.equals("C9")) {
                elementReceiverMap.put("C9", new AtomType(Carbon.INSTANCE, "C9"));
            }else if (element.equals("H8A")) {
                elementReceiverMap.put("H8A", new AtomType(Hydrogen.INSTANCE, "H8A"));
            }else if (element.equals("H6A")) {
                elementReceiverMap.put("H6A", new AtomType(Hydrogen.INSTANCE, "H6A"));
            }else if (element.equals("H8B")) {
                elementReceiverMap.put("H8B", new AtomType(Hydrogen.INSTANCE, "H8B"));
            }else if (element.equals("H6B")) {
                elementReceiverMap.put("H6B", new AtomType(Hydrogen.INSTANCE, "H6B"));
            }else if (element.equals("Cu1")) {
                elementReceiverMap.put("Cu1", new AtomType(Copper.INSTANCE, "Cu1"));
            }else if (element.equals("H2")) {
                elementReceiverMap.put("H2", new AtomType(Hydrogen.INSTANCE, "H2"));
            }else if (element.equals("H4")) {
                elementReceiverMap.put("H4", new AtomType(Hydrogen.INSTANCE, "H4"));
            }else if (element.equals("H9A")) {
                elementReceiverMap.put("H9A", new AtomType(Hydrogen.INSTANCE, "H9A"));
            }else if (element.equals("H9B")) {
                elementReceiverMap.put("H9B", new AtomType(Hydrogen.INSTANCE, "H9B"));
            }else if (element.equals("H7A")) {
                elementReceiverMap.put("H7A", new AtomType(Hydrogen.INSTANCE, "H7A"));
            }else if (element.equals("H7B")) {
                elementReceiverMap.put("H7B", new AtomType(Hydrogen.INSTANCE, "H7B"));
            }else if (element.equals("H9C")) {
                elementReceiverMap.put("H9C", new AtomType(Hydrogen.INSTANCE, "H9C"));
            }
           /* switch (element){

                case "Cu":
                    elementReceiverMap.put("Cu", new AtomType(Copper.INSTANCE, "Cu"));
                case "H":
                    elementReceiverMap.put("H", new AtomType(Hydrogen.INSTANCE, "H"));
                case "C":
                    elementReceiverMap.put("C", new AtomType(Carbon.INSTANCE, "C"));
            }*/

        }
        return elementReceiverMap.get(element);
    }

    public ArrayList<String[]> parseLineReader(String filePath, double[] boxSize, double[] boxAngle, Map<String, Vector> positionCoordinatesMap,Map<String, String> atomTypeMap  ){

        ArrayList<String[]> symmetryOperations = new ArrayList<>();
        List<String> atomSiteColumns = new ArrayList<>();
        List<String[]> atomSiteData = new ArrayList<>();
        boolean inAtomSiteLoop = false;
        boolean captureNextLines = false;
        Pattern patternSize = Pattern.compile("_cell_length_[abc]\\s+([\\d.]+)\\(\\d+\\)");
        Pattern patternAngle = Pattern.compile("_cell_angle_(alpha|beta|gamma)\\s+(\\d+)");
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {

                Matcher matcherSize = patternSize.matcher(line);
                Matcher matcherAngle = patternAngle.matcher(line);
                if (line.trim().startsWith("loop_")) {
                    captureNextLines = false; // Reset on every new loop
                }
                if (line.trim().equals("_space_group_symop_operation_xyz")) {
                    captureNextLines = true; // Set flag to start capturing the next lines
                    continue;
                }
                if (captureNextLines && line.trim().startsWith("'") && line.trim().endsWith("'")) {
                    String trimmedLine = line.trim().substring(1, line.trim().length() - 1); // Remove leading and trailing '
                    String[] operations = trimmedLine.split(",");
                    for (int i = 0; i < operations.length; i++) {
                        operations[i] = operations[i].trim(); // Trim each operation
                    }
                    symmetryOperations.add(operations);
                }
                if (matcherSize.find()) {
                    double value = Double.parseDouble(matcherSize.group(1));
                    if (line.contains("_cell_length_a")) boxSize[0] = value;
                    else if (line.contains("_cell_length_b")) boxSize[1] = value;
                    else if (line.contains("_cell_length_c")) boxSize[2] = value;
                } else if (matcherAngle.find()) {
                    int value = Integer.parseInt(matcherAngle.group(2));
                    if (line.contains("_cell_angle_alpha")) boxAngle[0] = value;
                    else if (line.contains("_cell_angle_beta")) boxAngle[1] = value;
                    else if (line.contains("_cell_angle_gamma")) boxAngle[2] = value;
                }

                if (line.trim().equals("loop_")) {
                    inAtomSiteLoop = false; // Reset flag for each new loop
                    atomSiteColumns.clear(); // Clear previous column headers
                } else if (line.trim().startsWith("_atom_site_")) {
                    inAtomSiteLoop = true;
                    atomSiteColumns.add(line.trim());
                } else if (inAtomSiteLoop && !line.trim().startsWith("_")) {
                    // Process the atom site data line
                    String[] dataParts = line.trim().split("\\s+");
                    if (dataParts.length >= atomSiteColumns.size()) {
                        // Extract relevant data based on the column headers
                        atomSiteData.add(dataParts);
                    }
                } else if (!line.trim().startsWith("_") && inAtomSiteLoop) {
                    // Exit the loop once non-relevant lines are detected
                    inAtomSiteLoop = false;
                }

            /*   if (line.trim().equals("loop_")) {
                    inAtomSiteLoop = false; // Reset flag for new loop
                    readyToCaptureData = false; // Reset ready flag
                } else if (line.trim().startsWith("_atom_site_label")) {
                    inAtomSiteLoop = true; // Start of atom site data loop detected
                } else if (inAtomSiteLoop && (line.trim().equals("_atom_site_disorder_group") || line.trim().equals("_atom_site_disorder_assembly"))) {
                    // After reading the specific headers, set ready to capture data
                    readyToCaptureData = true;
                } else if (readyToCaptureData && !line.trim().startsWith("_")) {
                    // Process the atom site data line
                    String[] dataParts = line.trim().split("\\s+");
                    if (dataParts.length > 4) { // Ensure there are enough parts for coordinates
                        atomSiteData.add(dataParts);
                    }
                }*/
            }
            System.out.println("Done");
        } catch (IOException e) {
            e.printStackTrace();
        }

        for (String[] atomData : atomSiteData) {
            System.out.println(Arrays.toString(atomData));
            // Assuming the first three positions after the label and type are x, y, z coordinates
            // and the atom label is always the first data piece.
            String atomLabel = atomData[0];
            String atomType = atomData[1];
            atomTypeMap.put(atomLabel, atomType);

            // Extract coordinates before the "("
            double x = Double.parseDouble(atomData[2].split("\\(")[0]);
            double y = Double.parseDouble(atomData[3].split("\\(")[0]);
            double z = Double.parseDouble(atomData[4].split("\\(")[0]);
            Vector positn = Vector.of(x, y, z);
            positionCoordinatesMap.put(atomLabel, positn);

            System.out.println(atomLabel + " " + atomType + " (" + x + ", " + y + ", " + z + ")");
        }

        return symmetryOperations;
    }


    public static ArrayList<double[]> getCoordinates(double[] coordinatesMain, ArrayList<String[]> symmentryOperation , double[] boxSize, double multiplier){
        ArrayList<double[]> positionCoordinates = new ArrayList<>();
        for (int j=0; j<symmentryOperation.size(); j++){
            double[] replicaCoordinates = new double[3];
            String[] symmetryOperationArray = symmentryOperation.get(j);
            for (int i = 0; i < symmetryOperationArray.length; i++) {
                //   System.out.println(Arrays.toString(symmetryOperationArray));
                switch (symmetryOperationArray[i]) {
                    case "x":
                        replicaCoordinates[i] = multiplier*coordinatesMain[0];
                        break;
                    case "y":
                        replicaCoordinates[i] =multiplier*coordinatesMain[1];
                        break;
                    case "z":
                        replicaCoordinates[i] = multiplier*coordinatesMain[2];
                        break;
                    case "-x":
                        replicaCoordinates[i] = -multiplier*coordinatesMain[0];
                        break;
                    case "-y":
                        replicaCoordinates[i] =-multiplier*coordinatesMain[1];
                        break;
                    case "-z":
                        replicaCoordinates[i] = -multiplier*coordinatesMain[2];
                        break;
                    case "x+1/2":
                        replicaCoordinates[i] =multiplier*(coordinatesMain[0]+0.5*boxSize[0]) ;
                        break;
                    case "y+1/2":
                        replicaCoordinates[i] =multiplier * (coordinatesMain[1]+0.5*boxSize[1]);
                        break;
                    case "z+1/2":
                        replicaCoordinates[i] = multiplier*(coordinatesMain[2]+0.5*boxSize[2]);
                        break;
                    case "x-1/2":
                        replicaCoordinates[i] = multiplier*(coordinatesMain[0]-0.5*boxSize[0]);
                        break;
                    case "y-1/2":
                        replicaCoordinates[i] =multiplier*(coordinatesMain[1]-0.5*boxSize[1]);
                        break;
                    case "z-1/2":
                        replicaCoordinates[i] = multiplier*(coordinatesMain[2]-0.5*boxSize[2]);
                        break;
                    case "-x-1/2":
                        replicaCoordinates[i] = multiplier*(-coordinatesMain[0]-0.5*boxSize[0]);
                        break;
                    case "-y-1/2":
                        replicaCoordinates[i] =multiplier*(-coordinatesMain[1]-0.5*boxSize[1]);
                        break;
                    case "-z-1/2":
                        replicaCoordinates[i] = multiplier*(-coordinatesMain[2]-0.5*boxSize[2]);
                        break;
                    case "-x+1/2":
                        replicaCoordinates[i] = multiplier*(-coordinatesMain[0]+0.5*boxSize[0]);
                        break;
                    case "-y+1/2":
                        replicaCoordinates[i] =multiplier*(-coordinatesMain[1]+0.5*boxSize[1]);
                        break;
                    case "-z+1/2":
                        replicaCoordinates[i] = multiplier*(-coordinatesMain[2]+0.5*boxSize[2]);
                        break;
                    /*default:
                        System.out.println("Not found" + symmetryOperationArray[i]);*/
                        // Add cases here for other transformation rules as needed
                }
               // System.out.println(Arrays.toString(replicaCoordinates));
            }
            positionCoordinates.add(replicaCoordinates);
        }
        return positionCoordinates;
    }

    public static void main(String[] args) {
        String confNameOne = "F://Avagadro//mopstrands//new//cif01";
        cifReader cifReader = new cifReader();
        cifReader.readCIFFile(confNameOne, 5);
        cifReader.speciesCIF(confNameOne, false);
        //System.exit(1);
        ArrayList<String[]> symmetryOperations = new ArrayList<>();
        Double[] boxSize = new Double[]{18.948, 18.948,18.948};
        Double[] coordinatesMain = new Double[]{1.0,2.0,3.0};
       /* symmetryOperations.add(new String[]{"x-1/2", "y+1/2", "z"});
        symmetryOperations.add(new String[]{"x+1/2", "y-1/2", "z-1/2"});
        symmetryOperations.add(new String[]{"x", "y", "z+1/2"});
        symmetryOperations.add(new String[]{"-y+1/2", "x", "z"});
        symmetryOperations.add(new String[]{"y", "-x+1/2", "z"});
        symmetryOperations.add(new String[]{"x+1/2", "-y", "-z"});
        symmetryOperations.add(new String[]{"-x", "y+1/2", "-z"});
        symmetryOperations.add(new String[]{"-x+1/2", "-y+1/2", "z"});
        symmetryOperations.add(new String[]{"y+1/2", "x+1/2", "-z"});
        symmetryOperations.add(new String[]{"-y", "-x", "-z"});
        symmetryOperations.add(new String[]{"-x", "-y", "-z"});
        symmetryOperations.add(new String[]{"y-1/2", "-x", "-z"});
        symmetryOperations.add(new String[]{"-y", "x-1/2", "-z"});
        symmetryOperations.add(new String[]{"-x-1/2", "y", "z"});
        symmetryOperations.add(new String[]{"x", "-y-1/2", "z"});
        symmetryOperations.add(new String[]{"x-1/2", "y-1/2", "-z"});
        symmetryOperations.add(new String[]{"-y-1/2", "-x-1/2", "z"});
        symmetryOperations.add(new String[]{"x", "-y", "z"});
        symmetryOperations.add(new String[]{"y", "x", "z"});
       //getCoordinates("C", positionList, coordinatesMain, symmetryOperations, boxSize );
       /* for (int i=0; i<symmetryOperations.size(); i++){
            Double[] coordinates = getCoordinates(coordinatesMain, symmetryOperations.get(i), boxSize);
          //  System.out.println(Arrays.toString(coordinates));
        }*/

        //System.out.println(Arrays.deepToString(symmetryOperations.toArray()));
    }
}

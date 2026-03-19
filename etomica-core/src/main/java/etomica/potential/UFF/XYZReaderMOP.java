package etomica.potential.UFF;

import etomica.atom.AtomType;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class XYZReaderMOP {
    public ISpecies speciesMOP;
    public Map<Integer, String> atomMap = new HashMap<>();
    public Map<Integer, AtomType> atomTypeMap = new HashMap<>();
    public Map<String, AtomType> typeMap = new HashMap<>();
    public Map<Integer, Vector> positions = new HashMap<>();
    public ISpecies getSpeciesMOP(String confName, boolean isInfinite, boolean isDynamic){
        XYZReaderMOP xyzReaderMOP = new XYZReaderMOP();
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        SpeciesBuilder speciesBuilderNew =  new SpeciesBuilder(Space3D.getInstance());
        pdbReaderMOP.readPDBFileCrystal(confName);
        atomMap = pdbReaderMOP.getAtomMap();
        positions = pdbReaderMOP.getPositions();
        Vector centreMOP = new Vector3D(0.0, 0.0, 0.0);
        AtomType typeNew;
        Vector dr = Vector.d(centreMOP.getD());
        Vector position = new Vector3D();
        for (int i = 0; i < atomMap.size(); i++){
            String symbol = atomMap.get(i);
            AtomType atomName = pdbReaderMOP.returnElement(symbol, isInfinite);
            atomTypeMap.put(i, atomName);
            if (typeMap.containsKey(symbol)) {
                typeNew = typeMap.get(symbol);
            } else {
                typeNew = atomName;
                typeMap.put(symbol, typeNew);
            }
            position = positions.get(i);
            speciesBuilderNew.addAtom(typeNew, position,  "");
        }
       // System.out.println(typeMap);
       // System.out.println(atomMap);
        speciesMOP= speciesBuilderNew.setDynamic(isDynamic).build();
        return speciesMOP;
    }
    public void readXYZFile(String confName ){
        String fileName = confName+".xyz";
        FileReader fileReader;
        Map<String, AtomType> typeMap = new HashMap<>();
        ArrayList<Integer> currentAtomList = null;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try{
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line ;
            Integer atomNumber = 0;
            while ((line = bufferedReader.readLine()) != null){
                parseLineReader(line, typeMap, atomMap, positions);
            }
            fileReader.close();


        }catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        setPositions(positions);
        setAtomMap(atomMap);
    }

    protected  void parseLineReader (String line, Map<String, AtomType> typeMap, Map<Integer, String> atomMap, Map<Integer, Vector> positionMap){
        line = line.trim();
        String[] parts = line.split(" ");
        if (parts.length == 4 && parts[0].length() < 2 && parts[1].length() > 10){
            double x = Double.parseDouble(parts[1]);
            double y = Double.parseDouble(parts[2]);
            double z = Double.parseDouble(parts[3]);
            Vector posn = Vector.of(x, y, z);
            String symbol = parts[0];
            int atomNumber = Integer.parseInt(parts[4]);
            positionMap.put(atomNumber, posn);
            atomMap.put(atomNumber, symbol);

        }
    }

    private void setPositions(Map<Integer, Vector> positions){this.positions = positions;}
    public Map<Integer, Vector> getPositions(){return positions;}
    private void setAtomMap(Map<Integer, String> atomMap){this.atomMap = atomMap;}
    public Map<Integer, String> getAtomMap(){return atomMap;}
}

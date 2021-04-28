package etomica.config;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.parser.ParserAMBER;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class ConformationFile {
    public static ConformationGeneric makeConformation(Space space, String fileName) {
        Vector[] coords = null;
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        List<String> lines = new ArrayList<String>();
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line = null;
            while ((line = bufReader.readLine()) != null) {
                lines.add(line.trim());
            }
            bufReader.close();
            }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
            //System.out.println(lines);
            //return makeStuffFromLines(lines, opts);
        int i = 0;
        Iterator<String> linesIterator = lines.iterator();
        while (linesIterator.hasNext()) {
            String a = linesIterator.next();
            String[] fields = a.split("\\s+");
            if (i == 0){
                coords = new Vector[Integer.parseInt(fields[0])];
                ++i;
            }
            else{
                if(!a.isEmpty()){
                    coords[i-1] = space.makeVector();
                    for (int j =1; j < fields.length; ++j) {
                        coords[i-1].setX(j - 1, Double.parseDouble(fields[j]));
                    }
                    ++i;
                }
            }
        }

        ConformationGeneric conformationGeneric = new ConformationGeneric(coords);
        return conformationGeneric;
        //ConformationFile.readConfiguration("CCCCCC.prmtop", new ParserAMBER.Options());
    }
}

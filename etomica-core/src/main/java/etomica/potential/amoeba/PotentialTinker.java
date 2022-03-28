/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.space.Vector;
import etomica.units.*;

import java.io.*;
import java.util.Map;

public class PotentialTinker {

    protected final Box box;
    public static String tinkerPath = "/home/andrew/build/tinker/";
    protected final String execPath, forceFieldPath;
    protected final int[][][] bonding;
    protected final Map<AtomType,Integer> typeMap;

    public PotentialTinker(Box box, String forceField, int[][][] bonding, Map<AtomType,Integer> typeMap) {
        this.box = box;
        execPath = tinkerPath+"source/analyze.x";
        forceFieldPath = tinkerPath+"params/"+forceField+".prm";
        this.bonding = bonding;
        this.typeMap = typeMap;
    }

    protected void writeXyz(IMolecule... molecules) {
        int numAtoms = 0;
        for (IMolecule m : molecules) {
            numAtoms += m.getChildList().size();
        }
        try {
            FileWriter fw = new FileWriter("molecules.xyz");

            fw.write(""+numAtoms+"\n");

            int ia = 1;
            for (IMolecule m : molecules) {
                int is = m.getType().getIndex();
                for (IAtom a : m.getChildList()) {
                    Vector p = a.getPosition();
                    fw.write(ia + " " + a.getType().getElement().getSymbol() + " " + p.getX(0) + " " + p.getX(1) + " " + p.getX(2) + " " + typeMap.get(a.getType()));
                    int[] b = bonding[is][a.getIndex()];
                    for (int j = 0; j < b.length; j++) {
                        fw.write(" " + b[j]);
                    }
                    fw.write("\n");
                }
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public double energy(IMolecule... molecules) {
        writeXyz(molecules);
        return runTinker();
    }

    protected double runTinker() {
        double u = Double.NaN;
        try{
            File file = new File("tinker.out");
            if (file.exists()) {
                if (!file.delete()) {
                    throw new RuntimeException("unable to clear out tinker.out");
                }
            }
            ProcessBuilder pb = new ProcessBuilder(execPath, "molecules.xyz");
            pb.redirectErrorStream(true);
            Process proc = pb.start();
            OutputStream out = proc.getOutputStream();
            out.write("E".getBytes());
            BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));

            String line;
            Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("Total Potential Energy :")) {
                    String[] bits = line.split("\\s+");
                    double ukcalpmole = Double.parseDouble(bits[4]);
                    u = kcalpmole.toSim(ukcalpmole);
                }
            }

            proc.waitFor();
        }
        catch (IOException | InterruptedException e){
            System.out.println("Problem running Tinker");
            throw new RuntimeException(e);
        }
        return u;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
}

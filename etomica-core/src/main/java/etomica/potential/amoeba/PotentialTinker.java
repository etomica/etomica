/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.space.Vector;
import etomica.units.*;
import etomica.util.collections.IntArrayList;

import java.io.*;
import java.util.Map;

public class PotentialTinker implements IPotentialMolecular {

    public static String tinkerPath = "/home/andrew/build/tinker0/";
    protected final String execPath, forceFieldPath;
    protected final IntArrayList[][] bonding;
    protected final Map<AtomType,Integer> typeMap;

    public PotentialTinker(String forceField, IntArrayList[][] bonding, Map<AtomType,Integer> typeMap) {
        execPath = tinkerPath+"source/analyze.x";
        forceFieldPath = tinkerPath+"params/"+forceField+".prm";
        this.bonding = bonding;
        this.typeMap = typeMap;
        writeKey();
    }

    protected void writeKey() {
        try {
            FileWriter fw = new FileWriter("molecules.key");
            fw.write("parameters "+forceFieldPath+"\n");
            fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
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
            int offset = 1;
            for (IMolecule m : molecules) {
                int is = m.getType().getIndex();
                for (IAtom a : m.getChildList()) {
                    Vector p = a.getPosition();
                    fw.write(ia + " " + a.getType().getElement().getSymbol() + " " + p.getX(0) + " " + p.getX(1) + " " + p.getX(2) + " " + typeMap.get(a.getType()));
                    IntArrayList b = bonding[is][a.getIndex()];
                    for (int j = 0; j < b.size(); j++) {
                        fw.write(" " + (b.getInt(j)+offset));
                    }
                    fw.write("\n");
                    ia++;
                }
                offset += m.getChildList().size();
            }
            fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public double energy(IMoleculeList molecules) {
        IMolecule[] array = new IMolecule[molecules.size()];
        for (int i=0; i<molecules.size(); i++) {
            array[i] = molecules.get(i);
        }
        return energy(array);
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
            ProcessBuilder pb = new ProcessBuilder(execPath, "molecules.xyz", "e");
            pb.redirectErrorStream(true);
            Process proc = pb.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));

            String line;
            Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
            while ((line = reader.readLine()) != null) {
                if (line.matches(".Total Potential Energy : .* Kcal/mole.*")) {
//                if (line.startsWith("Total Potential Energy :")) {
                    String[] bits = line.split("\\s+");
                    double ukcalpmole = Double.parseDouble(bits[5]);
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

    public static class P1Tinker extends PotentialTinker implements IPotentialMoleculeSingle {
        public P1Tinker(String forceField, IntArrayList[][] bonding, Map<AtomType, Integer> typeMap) {
            super(forceField, bonding, typeMap);
        }

        @Override
        public double energy(IMolecule molecule) {
            return super.energy(molecule);
        }
    }
}

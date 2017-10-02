/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.box.Box;
import etomica.config.Configuration;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ConfigurationLammps implements Configuration {
    protected final Space space;
    protected final String filename;
    protected final ISpecies[] species;
    protected Vector shift;
    protected double Lxy;
    protected double topPadding;
    protected double zMax;
    
    public ConfigurationLammps(Space space, String filename, ISpecies topWall, ISpecies bottomWall, ISpecies fluid) {
        this.space = space;
        this.filename = filename;
        this.species = new ISpecies[]{topWall, fluid, bottomWall};
        shift = space.makeVector();
    }
    
    public void setTopPadding(double newTopPadding) {
        topPadding = newTopPadding;
    }
    
    public void initializeCoordinates(Box box) {
        List<Vector>[] coords = new List[3];
        coords[0] = new ArrayList<Vector>();
        coords[1] = new ArrayList<Vector>();
        coords[2] = new ArrayList<Vector>();
        double zMin = Double.POSITIVE_INFINITY;
        zMax = -Double.POSITIVE_INFINITY;
        boolean readCoords = false;
        try {
            FileReader fr = new FileReader(filename);
            BufferedReader bufReader = new BufferedReader(fr);
            String line = null;
            while ((line = bufReader.readLine()) != null) {
                String[] bits = line.split("[\\t ]");
                if (bits.length == 4 && bits[2].equals("xlo")) {
                    Lxy = Double.parseDouble(bits[1]);
                    continue;
                }
                if (line.matches("Atoms # atomic")) {
                    readCoords = true;
                    continue;
                }
                
                if (!readCoords || bits.length < 5) continue;
                int aType = Integer.parseInt(bits[1]);
                Vector xyz = space.makeVector();
                for (int i=0; i<3; i++) {
                    xyz.setX(i, Double.parseDouble(bits[2+i]));
                }
                if (xyz.getX(2) > zMax) zMax = xyz.getX(2);
                if (xyz.getX(2) < zMin) zMin = xyz.getX(2);
                coords[aType-1].add(xyz);
            }
            bufReader.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        double Lz = zMax-zMin + 0.02 + topPadding;
        shift.setX(0, Lxy);
        shift.setX(1, Lxy);
        shift.setX(2, Lz);
        box.getBoundary().setBoxSize(shift);
        shift.TE(-0.5);
        shift.setX(2, -0.5*Lz - zMin);
        for (int i=0; i<3; i++) {
            box.setNMolecules(species[i], coords[i].size());
            IMoleculeList m = box.getMoleculeList(species[i]);
            for (int j=0; j<coords[i].size(); j++) {
                Vector p = m.getMolecule(j).getChildList().getAtom(0).getPosition();
                p.Ev1Pv2(coords[i].get(j), shift);
            }
        }
        zMax += shift.getX(2);
    }
    
    public Vector getShift() {
        return shift;
    }
    
    public double getLxy() {
        return Lxy;
    }
}

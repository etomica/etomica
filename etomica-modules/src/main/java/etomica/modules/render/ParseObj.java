/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.render;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import etomica.space3d.Vector3D;

public class ParseObj {
    
    public ParseObj(String fileName) {
        
        FileReader fileReader;
        
        int countV = 0;
        int bondedCount = 0;
        try {
            fileReader = new FileReader(fileName);
        } catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        
        double xmax = 0, xmin = 0, ymax = 0, ymin = 0;
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String currentLine;
//            for (int iLeaf=0; iLeaf<30; iLeaf++) {
            int[] indexRef = new int[100];
            double tol = 1e-14;
            while((currentLine = bufReader.readLine()) != null) {
                String[] coordStr = currentLine.split(" ");
                if(coordStr[0].equals("v") /*&& countV < 1000*/) {
                    Vector3D vertex = new Vector3D(Double.valueOf(coordStr[1]).doubleValue(),Double.valueOf(coordStr[2]).doubleValue(),Double.valueOf(coordStr[3]).doubleValue());
                    if(countV == indexRef.length) indexRef = expand(indexRef);
                    
                    //look for vertex position in list of ones read previously
                    boolean foundVertex = false;
                    for(int i=0; i<vertices.size(); i++) {
                        double r2 = vertex.Mv1Squared(vertices.get(i));
                        if(r2 < tol) {//found vertex in list
                            foundVertex = true;
                            indexRef[countV] = i;
                            break;
                        }
                    }
                    xmax = Math.max(xmax, vertex.getX(0));
                    xmin = Math.min(xmin, vertex.getX(0));
                    ymax = Math.max(ymax, vertex.getX(1));
                    ymin = Math.min(ymin, vertex.getX(1));
                    if(!foundVertex) {
                        indexRef[countV] = vertices.size();
                        vertices.add(vertex);
                    }
                    countV++;
                }
                if(coordStr[0].equals("vn")) {
                    //vertex normal; do not use
                }
                if(coordStr[0].equals("f")) {
                    int i0 = indexRef[Integer.valueOf(coordStr[1].split("//")[0]).intValue() - 1];
                    int i1 = indexRef[Integer.valueOf(coordStr[2].split("//")[0]).intValue() - 1];
                    int i2 = indexRef[Integer.valueOf(coordStr[3].split("//")[0]).intValue() - 1];
                    Vector3D v0 = vertices.get(i0);                    
                    Vector3D v1 = vertices.get(i1);                    
                    Vector3D v2 = vertices.get(i2);
                    if(!isBonded(i0,i1)) {
                        bondList.add(new BondInfo(i0, i1, v0.Mv1Squared(v1)));
                    } else bondedCount++;
                    if(!isBonded(i0,i2)) {
                        bondList.add(new BondInfo(i0, i2, v0.Mv1Squared(v2)));
                    } else bondedCount++;

                    if(!isBonded(i1,i2)) {
                        bondList.add(new BondInfo(i1, i2, v1.Mv1Squared(v2)));
                    } else bondedCount++;
                }
            }
            System.out.println("xmax, xmin: "+xmax+" "+xmin);
            System.out.println("ymax, ymin: "+ymax+" "+ymin);
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        
        doDeleteP(0,0.93*xmax);
        doDeleteM(0,0.93*xmin);
        doDeleteP(1,0.93*ymax);
        doDeleteM(1,0.93*ymin);
        nAtoms = vertices.size();
        System.out.println(nAtoms+" atoms from "+countV+" vertices; ");
        System.out.println(bondList.size()+" bonds after eliminating "+bondedCount+" redundant bonds");
    }
    
    private void doDeleteP(int j, double x) {
        boolean deleting = true;
        while(deleting) {
            deleting = false;
            for(int i=0; i<vertices.size(); i++) {
                if(vertices.get(i).getX(j) >= x) {
                    System.out.println("Deleting "+vertices.get(i).getX(0));
                    delete(i);
                    deleting = true;
                }
            }
        }
    }
    
    private void doDeleteM(int j, double x) {
        boolean deleting = true;
        while(deleting) {
            deleting = false;
            for(int i=0; i<vertices.size(); i++) {
                if(vertices.get(i).getX(j) <= x) {
                    System.out.println("Deleting "+vertices.get(i).getX(0));
                    delete(i);
                    deleting = true;
                }
            }
        }
    }

    
    /**
     * Removes the vertex referenced by the given index from the list of vertices,
     * removes its bonds from the list of bonds, shifts all vertex indexes down
     * to fill its slot, and returns a list of indexes (after shifting) of vertices
     * that were bonded to the given vertex.
     */
    protected ArrayList<Integer> delete(int k) {
        vertices.remove(k);

        //identify the vertex's bonds
        ArrayList<BondInfo> bonds = neighborBonds(k);
        
        //shift vertex indices above (and not including) deleted vertex
        int nBonds = bondList.size();
        for(int i=0; i<nBonds; i++) {
            bondList.get(i).shift(k);
        }

        //delete bonds and take note of bonding partners (with the new index) 
        ArrayList<Integer> nbrs = new ArrayList<Integer>(); 
        for(int i=0; i<bonds.size(); i++) {
            BondInfo bond = bonds.get(i);
            bondList.remove(bond);
            if(bond.i0 == k) nbrs.add(new Integer(bond.i1));
            else nbrs.add(new Integer(bond.i0));
        }
        return nbrs;
    }
    
    protected ArrayList<Integer> neighbors(int k) {
        ArrayList<Integer> neighbors = new ArrayList<Integer>();
        int nBonds = bondList.size();
        for(int i=0; i<nBonds; i++) {
            BondInfo bond = bondList.get(i);
            if(bond.i0 == k) neighbors.add(new Integer(bond.i1));
            if(bond.i1 == k) neighbors.add(new Integer(bond.i0));
        }
        return neighbors;
    }

    protected ArrayList<BondInfo> neighborBonds(int k) {
        ArrayList<BondInfo> bonds = new ArrayList<BondInfo>();
        int nBonds = bondList.size();
        for(int i=0; i<nBonds; i++) {
            BondInfo bond = bondList.get(i);
            if(bond.i0 == k || bond.i1 == k) bonds.add(bond);
        }
        return bonds;
    }
    
    private boolean isBonded(int i0, int i1) {
        for(int i=0; i<bondList.size(); i++) {
            BondInfo bond = bondList.get(i);
            if((i0 == bond.i0 && i1 == bond.i1) || (i1 == bond.i0 && i0 == bond.i1)) {
                return true;
            }
        }
        return false;
    }
    
    private int[] expand(int[] index) {
        int[] newIndex = new int[index.length+100];
        System.arraycopy(index, 0, newIndex, 0, index.length);
        return newIndex;
    }
    
    public int nAtoms;
    public ArrayList<BondInfo> bondList = new ArrayList<BondInfo>();
    ArrayList<Vector3D> vertices = new ArrayList<Vector3D>();
    
    public class BondInfo {
        public BondInfo(int i0, int i1, double bondLengthSquared) {
            this.i0 = i0;
            this.i1 = i1;
            this.bondLengthSquared = bondLengthSquared;
        }
        protected void shift(int i) {
            if(i0 > i) i0--;
            if(i1 > i) i1--;
        }
        int i0, i1;
        double bondLengthSquared;
    }
    public static void main(String[] args) {
        String file = "/Users/kofke/Documents/workspace/car self-assembly/mustang.txt";
        ParseObj parser = new ParseObj(file);

    }

}

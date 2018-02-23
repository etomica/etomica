/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Kelvin;
import etomica.util.Constants;

import java.io.*;
import java.net.URL;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

public class P3HydrogenManzhos implements IPotential{
    public static void main(String[] args) {       

        Space space = Space3D.getInstance();
        P3HydrogenManzhos pot1 = new P3HydrogenManzhos(space);
        double E = 0.0;
        Vector[] v2 = new Vector[6];
        int maxi = 10000;
        boolean garb = true;        
        try {
            FileWriter file = new FileWriter(garb?"GarberoglioM.dat":"ManzhosBCM.dat");                            
            for (int i=0; i<maxi; i++) {
                double rMove = 12.00*i/maxi;
                double rh2h2 = 0.0;
                if (garb) {
                    for (int j=0; j<6; j++) {
                        v2[j] = space.makeVector();
                        v2[j].E(0);
                    }
                    Vector hbl = space.makeVector();
                    hbl.Ev1Mv2(eqPos[0], eqPos[1]);
                    double hb = Math.sqrt(hbl.squared())/2.0;
                    hbl.E(0);
                    hbl.setX(2, hb);
                    v2[0].PE(hbl);
                    hbl.TE(-2);
                    v2[1].Ev1Pv2(v2[0], hbl);
                    v2[2].setX(0, rMove);
                    hbl.TE(-0.5);
                    v2[2].PE(hbl);
                    hbl.TE(-2);
                    v2[3].Ev1Pv2(v2[2], hbl);
                    v2[4].setX(0, 1);
                    v2[4].setX(1, Math.sqrt(3));
                    v2[4].normalize();
                    v2[4].TE(rMove);
                    hbl.TE(-0.5);
                    v2[4].PE(hbl);
                    hbl.TE(-2);
                    v2[5].Ev1Pv2(v2[4], hbl);
                    rh2h2 = rMove;

                }
                else {
                    for (int j=0; j<6; j++) {
                        v2[j] = space.makeVector();
                        v2[j].E(eqPos[j]);
                    }
                    Vector com0 = space.makeVector();
                    Vector com1 = space.makeVector();
                    Vector com2 = space.makeVector();
                    Vector n0 = space.makeVector();
                    Vector hbl = space.makeVector();
                    hbl.Ev1Mv2(v2[4], v2[5]);
                    double hb = Math.sqrt(hbl.squared())/2.0;
                    hbl.normalize();
                    hbl.TE(hb);
                    com0.Ev1Pv2(v2[0], v2[1]);
                    com0.TE(0.5);
                    com1.Ev1Pv2(v2[2], v2[3]);
                    com1.TE(0.5);
                    com2.Ev1Pv2(v2[4], v2[5]);
                    com2.TE(0.5);
                    com1.PE(com0);
                    com1.TE(0.5);
                    n0.Ev1Mv2(com2, com1);
                    n0.normalize();
                    n0.TE(rMove);
                    com2.Ev1Pv2(com1, n0);
                    v2[4].Ev1Pv2(com2, hbl);
                    hbl.TE(-2);
                    v2[5].Ev1Pv2(v2[4], hbl);
                    hbl.Ev1Mv2(com0, com2);
                    rh2h2 = Math.sqrt(hbl.squared());                    
                }
                E = pot1.vH2H2H2(v2);
                file.write((garb?rh2h2:BohrRadius.UNIT.fromSim(rh2h2))+" "+(garb?Kelvin.UNIT.fromSim((E*Constants.PLANCK_H*1E-8*Constants.LIGHT_SPEED)):E)+"\n");
            }                
            file.close();
        } catch(IOException e) {
            throw new RuntimeException("caught IOException: " + e.getMessage());            
        }
    }
    protected Boundary boundary;
    protected final static int D = 12;
    protected final static int d = 9;
    protected final static int L = 1;
    protected final static int N = 80;
    protected double[] y = new double[d];
    protected double[][] A = new double[d][D];
    protected double[][] w = new double[N][d];
    protected double[] cn = new double[N];
    protected double[] b = new double[d];
    protected double[] dn = new double[N];
    protected double[] xMinp = new double[D];
    protected double[] xMaxp = new double[D];
    protected double xMint = -30407.289000;
    protected double xMaxt = 2854.092000;
    protected final double d0 = 0.126575;
    protected final double lambda = BohrRadius.UNIT.toSim(4.00);
    protected static Vector[] eqPos = new Vector[6];
    protected Vector vec;
    protected double[][] xPos = {{0.0374,-0.2422,0.2792},{-0.0374,0.2422,-0.2792},{-0.016,-1.2877,2.9799},{-0.0196,-2.0073,2.7949},{-0.1047,1.5911,2.54},{0.063,1.7936,3.2349}};
    protected int[][] nPerm = {{1,2,3,4,5,6,7,8,9,10,11,12},{1,2,3,4,6,5,8,7,10,9,12,11},{2,1,4,3,5,6,7,8,11,12,9,10},{3,4,1,2,7,8,5,6,9,10,11,12},{5,6,7,8,1,2,3,4,9,11,10,12},{1,3,2,4,9,10,11,12,5,6,7,8},{9,11,10,12,5,7,6,8,1,3,2,4}};
    public P3HydrogenManzhos(Space space) {
        for (int i=0; i<6; i++) {            
            eqPos[i] = space.makeVector();
            eqPos[i].E(xPos[i]);
        }
        vec = space.makeVector();
        getData();        
    }
    public double vH2H2H2 (Vector[] v1) { // this method is for debugging only, the energy method doesn't call this
        double[] q0 = new double [D];
        double[] q = new double [D];
        double[] qScaled = new double[D];
        int k = 0;
        for (int i=0; i<4; i++) {
            for (int j=(i%2 == 0 ? i+2:i+1); j<6; j++) {
                vec.Ev1Mv2(v1[j], v1[i]);
                double Rij = Math.sqrt(vec.squared());
                q0[k] = Math.exp(-Rij/lambda);
                k++;
            }
        }
        for (int i=2; i<6; i++) {
            q[i] = q0[i];
        }
        q0[2] = q[4];
        q0[3] = q[5];
        q0[4] = q[2];
        q0[5] = q[3];
        double f = 0.00;
        int maxPerm = 7;
        for (int perm=0; perm<maxPerm; perm++) {
            for (int j=0; j<D; j++) {
                q[j] = q0[nPerm[perm][j]-1];
            }            
            // scaling input coordinates
            for (int i=0; i<D; i++) {
                qScaled[i] = -1 + 2*(q[i] - xMinp[i])/(xMaxp[i] - xMinp[i]);
            }

            for (int i=0; i<d; i++) {
                y[i] = b[i];
                for (int j=0; j<D; j++) {
                    y[i] += A[i][j]*qScaled[j];
                }
            }

            double ePerm = d0;
            for (int i=0; i<N; i++) {
                double x = dn[i];
                for (int j=0; j<d; j++) {
                    x += w[i][j]*y[j];
                }            
                ePerm += cn[i]*sigma(x);
            }

            ePerm = sigma(ePerm);
            ePerm = xMint + 0.5*(ePerm+1.0)*(xMaxt - xMint);
            f += ePerm/(1000.0*maxPerm);
        }
        return (f*1E-8*Constants.PLANCK_H*Constants.LIGHT_SPEED);
    }

    protected double sigma (double x) {
        return (-1 + 2/(1+Math.exp(-2*x)));               
    }
    protected void getData() {
        int t = 0;
        InputStreamReader isr = null;
        FileReader fileReader = null;
        String fileName = "P3HydrogenManzhos_allData.dat";
        URL url = getClass().getResource("/etomica/potential/" + fileName);
        try {            
            if (url.getProtocol().equals("jar")) {
                URL myURL = new URL(url.getPath());
                // myURL.getPath() will be something like
                // jar:file:/usr/users/bob/ex_vir_hs.jar!/etomica/graph/model/impl/graph6fmt.zip
                String jarFileName = myURL.getPath().split("!")[0];
                JarFile jarFile = new JarFile(jarFileName);
                // strip off leading slash
                JarEntry jarEntry = jarFile.getJarEntry(myURL.getPath().split("!")[1].substring(1));
                InputStream jis = jarFile.getInputStream(jarEntry);
                isr = new InputStreamReader(jis);
                t = 1;
            }
            else {              
                fileReader = new FileReader(fileName);
            }
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader((t == 0 || isr != null)?fileReader:isr);
            String[] str0 = bufReader.readLine().trim().split(" +");
            for (int i=0; i<D; i++) {                
                xMinp[i] = Double.valueOf(str0[i]).doubleValue();
            }
            String[] str1 = bufReader.readLine().trim().split(" +");
            for (int i=0; i<D; i++) {                
                xMaxp[i] = Double.valueOf(str1[i]).doubleValue();
            }
            for (int i=0; i<d; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                for (int j=0; j<D; j++) {
                    A[i][j] = Double.valueOf(str[j]).doubleValue();
                }
            }
            for (int i=0; i<N; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                for (int j=0; j<d; j++) {
                    w[i][j] = Double.valueOf(str[j]).doubleValue();
                }
            }
            for (int i=0; i<N; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                cn[i] = Double.valueOf(str[0]).doubleValue();
            }
            for (int i=0; i<d; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                b[i] = Double.valueOf(str[0]).doubleValue();
            }
            for (int i=0; i<N; i++) {
                String[] str = bufReader.readLine().trim().split(" +");
                dn[i] = Double.valueOf(str[0]).doubleValue();
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
    }


    public double getRange() {        
        return Double.POSITIVE_INFINITY;
    }


    public void setBox(Box box) {
        boundary = box.getBoundary();        
    }


    public int nBody() {
        return 3;
    }
    public static class P3HydrogenManzhosMolecular extends P3HydrogenManzhos implements IPotentialMolecular {
        protected Vector[] v = new Vector[6];
        public P3HydrogenManzhosMolecular(Space space) {
            super(space);     
            for (int i=0; i<6; i++) {
                v[i] = space.makeVector();
            }
        }

        public double energy(IMoleculeList molecules) {
            for (int i=0; i<molecules.getMoleculeCount(); i++) {                
                v[2*i] = molecules.getMolecule(i).getChildList().get(0).getPosition();
                v[2*i+1] = molecules.getMolecule(i).getChildList().get(1).getPosition();
            }
            double[] q0 = new double [D];
            double[] q = new double [D];
            double[] qScaled = new double[D];
            int k = 0;
            for (int i=0; i<4; i++) {
                for (int j=i+2; j<6; j++) {
                    vec.Ev1Mv2(v[i], v[j]);
                    double Rij = Math.sqrt(vec.squared());
                    q0[k] = Math.exp(-Rij/lambda);
                    k++;
                }
            }
            for (int i=2; i<6; i++) { //reordering the array to match up with manzhos' input
                q[i] = q0[i];
            }
            q0[2] = q[4];
            q0[3] = q[5];
            q0[4] = q[2];
            q0[5] = q[3];
            double f = 0.0;
            int maxPerm = 7;
            for (int perm=0; perm<maxPerm; perm++) {
                for (int j=0; j<D; j++) {
                    q[j] = q0[nPerm[perm][j]-1];
                }            
                //            scaling input coordinates
                for (int i=0; i<D; i++) {
                    qScaled[i] = -1 + 2*(q[i] - xMinp[i])/(xMaxp[i] - xMinp[i]);
                }

                for (int i=0; i<d; i++) {
                    y[i] = b[i];
                    for (int j=0; j<D; j++) {
                        y[i] += A[i][j]*qScaled[j];
                    }
                }
                double ePerm = d0;
                for (int i=0; i<N; i++) {
                    double x = dn[i];
                    for (int j=0; j<d; j++) {
                        x += w[i][j]*y[j];
                    }            
                    ePerm += cn[i]*sigma(x);
                }
                ePerm = sigma(ePerm);
                ePerm = xMint + 0.5*(ePerm+1.0)*(xMaxt - xMint);
                f += ePerm/(1000.0*maxPerm);
            }
            return (f*1E-8*Constants.PLANCK_H*Constants.LIGHT_SPEED);
        } 
    }    

}

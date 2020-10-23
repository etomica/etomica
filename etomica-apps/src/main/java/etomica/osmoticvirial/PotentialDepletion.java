package etomica.osmoticvirial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.Potential2;
import etomica.potential.Potential2Spherical;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 */
public class PotentialDepletion extends Potential2 implements Potential2Spherical {

    protected final Vector dr;
    protected Boundary boundary;
    protected List<double[]> rwList;
    protected FileReader fileReader;
    protected String fileName;

    public PotentialDepletion(Space space, String file) {

        super(space);
        fileName = file;
        rwList = new ArrayList<>();
        dr = space.makeVector();

        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }

        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String data;
            while((data = bufReader.readLine()) != null) {
                String[] splitData = data.split("[ \t]+");
                double[] rW = new double[2];
                rW[0] = Double.parseDouble(splitData[0]);
                rW[1] = Double.parseDouble(splitData[1]);
                rwList.add(rW);
            }

            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
     }

    @Override
    public double u(double r2) {
        double r = Math.sqrt(r2);
        if(r>2.5) return 0;
        if(r<1.02) return Double.POSITIVE_INFINITY;
        double[][] rInterpolation = new double[2][2];
        int lookup = (int)((r-1.02)/0.04);
        rInterpolation[0] = rwList.get(lookup);
        rInterpolation[1] = rwList.get(lookup+1);
        return (rInterpolation[1][1]-rInterpolation[0][1]) / (rInterpolation[0][0] - rInterpolation[1][0]) * (r - rInterpolation[1][0]) + r;
    }

    @Override
    public double energy(IAtomList pair) {
        IAtom atom0 = pair.get(0);
        IAtom atom1 = pair.get(1);

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }

    @Override
    public double getRange() {
        double[] range = rwList.get(rwList.size()-1);
        return range[0];
    }

    @Override
    public void setBox(Box box) {
        boundary = box.getBoundary();

    }
}


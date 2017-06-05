/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

/**
 * Main method to drive AkimaSplineSmoother
 * 
 * @author Andrew Schultz
 */
public class AkimaSplineSmootherMain {

    public static String usage() {
        return "usage: smoother [--help] [-d2fac f] [-d2dfac f] [-d3fac f] [-d3dfac f]\n" +
               "                [-maxiter m] [-start start.dat] [-o out.dat]\n" +
               "                [-bump bumpInterval] input.dat";
    }
    
    public static String help() {
        return "  d2fac, d2dfac, d3fac or d3dfac must be non-negative.  At least one must be\n"+
               "     positive.\n" +
               "  if start.dat is specified, it will be used as initial y values (instead of\n"+
               "     values from input.dat)\n"+
               "  if out.dat is specified, smoothed data will be written there, along with\n" +
               "     first and second derivatives to out.dat.dy and out.dat.dy2\n"+
               "  if bumpInterval is specified, then the various factors will be bumped up\n"+
               "     or down by 10 occasionally to try to find new minima.";
    }
    
    public static void main(String[] args) {
        String infile = "B7.dat.scaled";
        String startfile = "";
        String outfile = "";
        double d2fac = 0;
        double d2dfac = 0;
        double d3fac = 0;
        double d3dfac = 0;
        int maxiter = 100;
        int iarg = 0;
        int bumpInterval = 0;
        boolean gotInfile = false;
        try {
            while (iarg < args.length) {
                if (args[iarg].equals("--help")) {
                    System.out.println(usage());
                    System.out.println(help());
                    System.exit(0);
                }
                if (args[iarg].equals("-d2fac")) {
                    if (iarg+1 == args.length) {
                        throw new RuntimeException(usage());
                    }
                    d2fac = Double.parseDouble(args[iarg+1]);
                    iarg+=2;
                }
                else if (args[iarg].equals("-d2dfac")) {
                    if (iarg+1 == args.length) {
                        throw new RuntimeException(usage());
                    }
                    d2dfac = Double.parseDouble(args[iarg+1]);
                    iarg+=2;
                }
                else if (args[iarg].equals("-d3fac")) {
                    d3fac = Double.parseDouble(args[iarg+1]);
                    iarg+=2;
                }
                else if (args[iarg].equals("-d3dfac")) {
                    d3fac = Double.parseDouble(args[iarg+1]);
                    iarg+=2;
                }
                else if (args[iarg].equals("-maxiter")) {
                    maxiter = Integer.parseInt(args[iarg+1]);
                    iarg+=2;
                }
                else if (args[iarg].equals("-start")) {
                    startfile = args[iarg+1];
                    iarg+=2;
                }
                else if (args[iarg].equals("-o")) {
                    outfile = args[iarg+1];
                    iarg+=2;
                }
                else if (args[iarg].equals("-bump")) {
                    bumpInterval = Integer.parseInt(args[iarg+1]);
                    iarg+=2;
                }
                else if (args[iarg].charAt(0) == '-') {
                    throw new RuntimeException("unknown arg '"+args[iarg]+"'");
                }
                else {
                    if (gotInfile) {
                        throw new RuntimeException("Multiple input files not allowed");
                    }
                    gotInfile = true;
                    infile = args[iarg];
                    iarg++;
                }
            }
        }
        catch (RuntimeException e) {
            System.out.println(usage());
            throw new RuntimeException(e);
        }
        if (d2fac * d2dfac * d3fac * d3dfac == 0 || d2fac < 0 || d2dfac < 0 || d3fac < 0 || d3dfac < 0) {
            System.out.println(usage());
            throw new RuntimeException("d2fac, d2dfac, d3fac or d3dfac must be non-negative.  At least one must be positive.");
        }
        IRandom rand = new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray());
        double[] x;
        double[] y;
        double[] dy;
        if (true) {
            double[][] readx = new double[3][0];
            readFile(infile, readx);
            x = readx[0];
            y = readx[1];
            dy = readx[2];
        }
        AkimaSplineSmoother fitter = new AkimaSplineSmoother(rand);
        fitter.setD2fac(d2fac);
        fitter.setD2dfac(d2dfac);
        fitter.setD3fac(d3fac);
        fitter.setD3dfac(d3dfac);
        fitter.setInputData(x, y, dy);

        double totSqErr = fitter.calcErr(0, x.length-1);
        System.out.println("i "+totSqErr+" "+fitter.sumSqDy+" "+fitter.sumSqD2+" "+fitter.sumSqD2D+" "+fitter.sumSqD3);

        if (!startfile.equals("")) {
            double[][] readx = new double[2][0];
            readFile(startfile, readx);
            fitter.setSmoothedY(readx[1]);
            totSqErr = fitter.calcErr(0, x.length-1);
            System.out.println("s "+totSqErr+" "+fitter.sumSqDy+" "+fitter.sumSqD2+" "+fitter.sumSqD2D+" "+fitter.sumSqD3);
        }
        if (maxiter == 0) {
            System.exit(0);
        }

        double oldErr = 0;
        double bumpFac = 0;
        double[] yold = new double[y.length];
        maxiter = ((maxiter + 4*bumpInterval - 1) / (4*bumpInterval)) * 4*bumpInterval;
        for (int i=0; i<maxiter; i++) {
            if (bumpInterval > 0 && i % bumpInterval == 0) {
                if (i/bumpInterval % 4 == 0) {
                    System.arraycopy(fitter.y, 0, yold, 0, yold.length);
                    oldErr = fitter.calcErr(0, x.length-1);
                    boolean up = rand.nextInt(2) == 0;
                    if (i==0) {
                        up = false;
                    }
                    System.out.println(up);
                    bumpFac = 10;
                    if (!up) {
                        bumpFac = 0.1;
                    }
                    if (d2fac != 0) {
                        d2fac *= bumpFac;
                        System.out.println("d2fac => "+d2fac);
                        fitter.setD2fac(d2fac);
                    }
                    if (d2dfac != 0) {
                        d2dfac *= bumpFac;
                        System.out.println("d2dfac => "+d2dfac);
                        fitter.setD2dfac(d2dfac);
                    }
                    if (d3fac != 0) {
                        d3fac *= bumpFac;
                        System.out.println("d3fac => "+d3fac);
                        fitter.setD3fac(d3fac);
                    }
                    if (d3dfac != 0) {
                        d3dfac *= bumpFac;
                        fitter.setD3dfac(d3dfac);
                        System.out.println("d3dfac => "+d3dfac);
                    }
                }
                else if (i/bumpInterval % 4 == 1) {
                    d2fac /= bumpFac;
                    fitter.setD2fac(d2fac);
                    d2dfac /= bumpFac;
                    fitter.setD2dfac(d2dfac);
                    d3fac /= bumpFac;
                    fitter.setD3fac(d3fac);
                    d3dfac /= bumpFac;
                    fitter.setD3dfac(d3dfac);
                }
                else if (i/bumpInterval % 4 == 3) {
                    double newErr = fitter.calcErr(0, x.length-1);
                    System.out.println("overall "+oldErr+" => "+newErr);
                    if (oldErr < newErr) {
                        // restore old state
                        System.out.println("(rejected)");
                        System.arraycopy(yold, 0, fitter.y, 0, yold.length);
                        System.out.println(" => "+fitter.calcErr(0, x.length-1));
                    }
                }
            }
            fitter.doStep();
            System.out.println((i+1)+" "+fitter.totSqErr+" "+fitter.sumSqDy+" "+fitter.sumSqD2+" "+fitter.sumSqD2D+" "+fitter.sumSqD3);
            if (!outfile.equals("")) {
                try {
                    FileWriter fw = new FileWriter(outfile);
                    for (int j=0; j<x.length; j++) {
                        fw.write(x[j]+" "+fitter.y[j]+"\n");
                    }
                    fw.close();
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }
                writeD(fitter, outfile, 10);
            }
        }
        if (outfile.equals("")) {
            for (int j=0; j<x.length; j++) {
                System.out.println(x[j]+" "+fitter.y[j]);
            }
        }
        else {
            try {
                FileWriter fw = new FileWriter(outfile);
                for (int j=0; j<x.length; j++) {
                    fw.write(x[j]+" "+fitter.y[j]+"\n");
                }
                fw.close();
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
            writeD(fitter, outfile, 10);
        }
        System.out.println(fitter.totSqErr+" "+fitter.sumSqDy+" "+d2fac*fitter.sumSqD2+" "+d2dfac*fitter.sumSqD2D+" "+d3fac*fitter.sumSqD3);
    }
    
    public static void readFile(String infile, double[][]x) {
        FileReader fileReader;
        try {
            fileReader = new FileReader(infile);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+infile+", caught IOException: " + e.getMessage());
        }
        ArrayList<Double>[] xLists = new ArrayList[x.length];
        for (int i=0; i<x.length; i++) {
            xLists[i] = new ArrayList<Double>();
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            while (true) {
                String line = bufReader.readLine();
                if (line == null) {
                    break;
                }
                String[] xstr = line.split(" +");
                for (int i=0; i<x.length; i++) {
                    xLists[i].add(Double.parseDouble(xstr[i]));
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading d.dat, caught IOException: " + e.getMessage());
        }

        for (int j=0; j<x.length; j++) {
            x[j] = new double[xLists[j].size()];
            for (int i=0; i<xLists[j].size(); i++) {
                x[j][i] = xLists[j].get(i);
            }
        }
    }

    public static void writeD(AkimaSplineSmoother fitter, String outbase, int nSubPoints) {
        double[][] dy12 = fitter.getDy12(nSubPoints);
        double[]x = fitter.x;

        try {
            FileWriter dyfw = new FileWriter(outbase+".dy");
            FileWriter dy2fw = new FileWriter(outbase+".dy2");

            int N = x.length;

            for (int i=0; i<N-1; i++) {
                for (int j=0; j<nSubPoints; j++) {
                    double dx = j*(x[i+1] - x[i])/nSubPoints;
                    double ix = x[i]+dx;
                    dyfw.write(ix+" "+dy12[0][i*nSubPoints+j]+"\n");
                    dy2fw.write(ix+" "+dy12[1][i*nSubPoints+j]+"\n");
                }
                if (i==N-2) {
                    dyfw.write(x[i+1]+" "+dy12[0][(i+1)*nSubPoints]+"\n");
                    dy2fw.write(x[i+1]+" "+dy12[1][(i+1)*nSubPoints]+"\n");
                }
            }

            dyfw.close();
            dy2fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}

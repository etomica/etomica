package etomica.virial.simulations;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.IntSet;
import etomica.virial.MeterVirialEBinMultiThreaded.MyData;
import etomica.virial.MeterVirialEBinMultiThreaded;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialSQWBinMultiRecalcW {


    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.nPoints = 6;
            params.maxRunNumber = 2;
            params.tRatio = 37.0925119792516;
        }
        
        final int nPoints = params.nPoints;
        final double tRatio = params.tRatio;
        if (tRatio <= 0) {
            throw new RuntimeException("must specify tRatio");
        }

        MeterVirialEBinMultiThreaded meter = new MeterVirialEBinMultiThreaded(null, null, null);
        meter.setTRatio(tRatio);
        String[] rawFiles = new String[params.maxRunNumber+1];
        for (int i=1; i<=rawFiles.length; i++) {
            rawFiles[i-1] = params.runName+nPoints+"_run"+i+"_raw.dat";
        }
        System.out.print("reading...");
        long t1 = System.currentTimeMillis();
        meter.readData(rawFiles, nPoints);
        long t2 = System.currentTimeMillis();
        System.out.println(String.format(" %6.3f",(t2-t1)*0.001));
        
        System.out.print("writing...");
        t1 = System.currentTimeMillis();
        meter.writeData(params.runName+nPoints+"_new_raw.dat");
        t2 = System.currentTimeMillis();
        System.out.println(String.format(" %6.3f",(t2-t1)*0.001));

        t1 = System.currentTimeMillis();
        meter.recomputeWeights(meter.getAllMyData(), meter.getTotalCount(), false, nPoints);
        t2 = System.currentTimeMillis();
        System.out.println(String.format("recomputed weights... %6.3f",(t2-t1)*0.001));

        System.out.print("writing weights...");
        t1 = System.currentTimeMillis();
        meter.writeWeights(params.runName+nPoints+"_new_weights.dat");
        t2 = System.currentTimeMillis();
        System.out.println(String.format(" %6.3f",(t2-t1)*0.001));

        Map<IntSet,MyData> allMyData = meter.getAllMyData();
        List<IntSet> pvs = new ArrayList<IntSet>();
        pvs.addAll(allMyData.keySet());
        Collections.sort(pvs);
        double[] sum = new double[1+nPoints*(nPoints-1)/2];
        double[] sumErrStdev = new double[sum.length];
        long steps = meter.getTotalCount();
        long totalSampleCount = 0;
        long totalNotScreenedCount = 0;
        int nSets = 0;
        FileWriter fw = null;
        try {
            fw = new FileWriter(params.runName+nPoints+"_new.dat");
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        double[] E0a = new double[sum.length];
        double[] E0a2 = new double[sum.length];
        for (IntSet pv : pvs) {
            MyData amd = allMyData.get(pv);
            long c = amd.unscreenedCount;

            nSets++;
            totalNotScreenedCount += c;
            long sc = amd.sampleCount;
            try {
                fw.write(Arrays.toString(pv.v)+" "+sc+"/"+c+"\n");
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }

            totalSampleCount += sc;

            for (int i=0; i<1+nPoints*(nPoints-1)/2; i++) {
                
                if (sc == 0) {
                    try {
                        fw.write("   0  0\n");
                    }
                    catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    continue;
                }


                double avg = amd.getAvg(i);
                double var = amd.getVar(i);
                double err = Math.sqrt(var/sc);

                try {
                    fw.write("   "+avg+"  "+err+"\n");
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }
                sum[i] += c*avg;

                E0a[i] += c*avg;
                E0a2[i] += c*avg*avg;
                if (var>0) {
                    sumErrStdev[i] += var/sc*c*c;
                }
            }
        }
        try {
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println(nSets+" sets");
        System.out.println();

        for  (int i=0; i<sum.length; i++) {
            double E0ave = E0a[i]/steps;
            double E0 = E0a2[i]/steps - E0ave*E0ave;
            sumErrStdev[i] /= steps;
            sum[i] /= steps;
            double finalErr = Math.sqrt((sumErrStdev[i] + E0)/steps);
            System.out.print(String.format("%2d average: %21.14e   error: %11.5e   # var frac: %5.3f\n", i, sum[i], finalErr, E0/(sumErrStdev[i] + E0)));
        }

        System.out.println();
        System.out.println("number time fraction: "+steps/(steps + totalSampleCount*tRatio));
        System.out.println("fraction not screened: "+((double)totalNotScreenedCount)/steps);
        System.out.println("fraction measured: "+((double)totalSampleCount)/totalNotScreenedCount);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 7;
        public String runName = "sqw";
        public int maxRunNumber = -1;
        public double tRatio = 0;
    }
    
}

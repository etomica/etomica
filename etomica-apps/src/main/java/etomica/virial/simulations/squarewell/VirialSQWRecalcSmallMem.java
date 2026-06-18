/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.squarewell;

import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.MeterVirialEBinMultiThreaded;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialSQWRecalcSmallMem {


    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.nPoints = 6;
            params.maxRunNumber = 4;
            params.tRatio = 37.0925119792516;
        }
        
        final int nPoints = params.nPoints;
        final double tRatio = params.tRatio;
        if (tRatio <= 0) {
            throw new RuntimeException("must specify tRatio");
        }

        MeterVirialEBinMultiThreaded meter = new MeterVirialEBinMultiThreaded(null, null, null, nPoints);
        MeterVirialEBinMultiThreaded.setTRatio(tRatio);
        String singleFile = null;
        if (params.maxRunNumber>1) {
            String[] rawFiles = new String[params.maxRunNumber];
            for (int i=1; i<=rawFiles.length; i++) {
                rawFiles[i-1] = params.runName+nPoints+"_run"+i+"_raw.dat";
            }
            System.out.print("reading and writing...");
            long t1 = System.currentTimeMillis();
            meter.readWriteData(rawFiles, params.runName+nPoints+"_new_raw.dat", nPoints);
            long t2 = System.currentTimeMillis();
            System.out.println(String.format(" %6.3f",(t2-t1)*0.001));
            singleFile = params.runName+nPoints+"_new_raw.dat";
        }
        else {
            singleFile = params.runName+nPoints+"_run1_raw.dat";
        }
        
        System.out.print("reading without covariance...");
        long t1 = System.currentTimeMillis();
        meter.readData(new String[]{singleFile}, nPoints, true);
        long t2 = System.currentTimeMillis();
        System.out.println(String.format(" %6.3f",(t2-t1)*0.001));
        
        t1 = System.currentTimeMillis();
        MeterVirialEBinMultiThreaded.recomputeWeights(meter.getAllMyData(), meter.getTotalCount(), true, nPoints);
        t2 = System.currentTimeMillis();
        System.out.println(String.format("recomputed weights... %6.3f",(t2-t1)*0.001));

        System.out.print("writing weights...");
        t1 = System.currentTimeMillis();
        meter.writeWeights(params.runName+nPoints+"_new_weights.dat");
        t2 = System.currentTimeMillis();
        System.out.println(String.format(" %6.3f",(t2-t1)*0.001));

        System.out.print("re-reading and processing...\n");
        t1 = System.currentTimeMillis();
        MeterVirialEBinMultiThreaded.readProcessData(singleFile, nPoints);
        t2 = System.currentTimeMillis();
        System.out.println(String.format("time: %6.3f",(t2-t1)*0.001));
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

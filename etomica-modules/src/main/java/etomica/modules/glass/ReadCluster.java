package etomica.modules.glass;

import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

public class ReadCluster {
    public static void main(String[] args) {
        ReadClusterParams params = new ReadClusterParams();
        if (args.length == 0) {
            params.nClusters = 10;
            params.filename = "sfac.bin";
            params.clusterFileIn = "cluster.bin";
            params.maxIterations = 10;
        } else {
            ParseArgs.doParseArgs(params, args);
        }
        long t1 = System.nanoTime();
        DataClusterer dataClusterer = new DataClusterer(params.nClusters, new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()));
        dataClusterer.setClusterNeighborDistance(params.neighborClusterDistance);
        DataClusterReader.readClusterFile(dataClusterer, params.filename, params.interval);
        long t2 = System.nanoTime();
        System.out.println("file read: " + (t2 - t1) / 1e9);
        int maxIterations = params.maxIterations;
        if (params.clusterFileIn != null) {
            dataClusterer.readClusterFile(params.clusterFileIn);
        } else if (maxIterations < 1) {
            maxIterations = 10;
        }
        if (maxIterations > 0) {
            dataClusterer.setMaxIterations(maxIterations);
            if (params.clusterFileOut != null) {
                dataClusterer.setClusterFileOut(params.clusterFileOut);
            }
            dataClusterer.findClusters();
        } else {
            dataClusterer.computePopulation();
        }
        long t3 = System.nanoTime();
        System.out.println("clusters assigned: " + (t3 - t2) / 1e9);
        dataClusterer.writeGraph("G.dot");
    }

    public static class ReadClusterParams extends ParameterBase {
        public String filename;
        public int nClusters = 100;
        public double neighborClusterDistance = 1.5;
        public int maxIterations = 0;
        public String clusterFileIn, clusterFileOut;
        public int interval = 1;
    }
}

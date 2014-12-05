/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.operations.ComponentSubst;
import etomica.graph.operations.ComponentSubst.ComponentSubstParameters;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.Parameters;
import etomica.graph.operations.Unary;

/**
 * This operation splits each e-bond component into a sum of permuted
 * (exchanged) nodes.  Nodes within an exchanged group are connected by exc
 * bonds.
 * 
 * @author Andrew Schultz
 */
public class ExchangeSplit implements Unary {

    protected byte[][][] labels = new byte[0][0][0];

    public Set<Graph> apply(Set<Graph> argument, Parameters params) {
        assert(params instanceof ExchangeSplitParameters);
        ExchangeSplitParameters myParams = (ExchangeSplitParameters)params;
        if (labels.length == 0) {
            labels = new byte[][][]{{{}},{{0}}};
        }
        ComponentSubst compSubst = new ComponentSubst();
        for (byte i=2; i<=myParams.maxExchangeNodes; i++) {
            Graph e = GraphFactory.createGraph(i);
            for (byte j=0; j<i*(i-1)/2; j++) {
                e.putEdge(j);
                e.getEdge(j).setColor(myParams.eBond);
            }
            Set<Graph> substSet = new HashSet<Graph>();

            int prevSize = labels[i-1].length;
            if (labels.length < i+1) {
                labels = (byte[][][])etomica.util.Arrays.addObject(labels, new byte[prevSize*i][i]);
            }
            List<List<Integer>> exchanges = new ArrayList<List<Integer>>();
            HashMap<String,Integer> exchangeIndex = new HashMap<String,Integer>();
            List<Integer> numExchanges = new ArrayList<Integer>();
            for (byte j=0; j<i; j++) {
                for (int k=0; k<prevSize; k++) {
                    byte[] myLabels = labels[i][j*prevSize+k];
                    // we build up the permuted labels from the labels with fewer nodes
                    for (byte l=0; l<i-1; l++) {
                        myLabels[l] = labels[i-1][k][l];
                        if (myLabels[l] == j) {
                            myLabels[i-1] = l;
                            myLabels[l] = (byte)(i-1);
                        }
                    }
                    myLabels[i-1] = j;
                    // myLabels are now our set of permuted labels.
                    // analyze this set
                    Set<Byte> done = new HashSet<Byte>();
                    List<Integer> thisExchange = new ArrayList<Integer>();
                    while (done.size() < i) {
                        // look for each exchange gruop
                        for (byte l=0; l<i; l++) {
                            if (done.contains(l)) continue;
                            int size = 1;
                            done.add(l);
                            int m = l;
                            // walk the exchange group and determine its size
                            while (myLabels[m] != l) {
                                done.add(myLabels[m]);
                                m = myLabels[m];
                                size++;
                            }
                            thisExchange.add(size);
                        }
                    }
                    // see if we've encountered this set of exchanges before
                    java.util.Collections.sort(thisExchange);
                    String exchangeStr = thisExchange.toString();
                    Integer thisIndex = exchangeIndex.get(exchangeStr);
                    if (thisIndex == null) {
                        // a new one
                        exchangeIndex.put(exchangeStr, exchanges.size());
                        exchanges.add(thisExchange);
                        numExchanges.add(1);
                    }
                    else {
                        // increment the count for this one
                        numExchanges.set(thisIndex, numExchanges.get(thisIndex)+1);
                    }
                }
            }

            // create a component for each exchange group, using the appropriate coefficient
            for (int j=0; j<exchanges.size(); j++) {
                List<Integer> exchange = exchanges.get(j);
                Graph s = e.copy();
                byte nodeId1 = 0;
                int sum = 0;
                for (Integer x : exchange) {
                    sum += x-1;
                    byte nodeId2 = (byte)(nodeId1 + x - 1);
                    for (byte iNode = nodeId1; iNode < nodeId2; iNode++) {
                        for (byte jNode = (byte)(iNode+1); jNode <= nodeId2; jNode++) {
                            s.getEdge(iNode,jNode).setColor(myParams.excBond);
                        }
                    }
                    nodeId1 += x;
                }
                int num = numExchanges.get(j);
                if (myParams.negativeExchange && sum%2==1) {
                    num = -num;
                }
                s.coefficient().multiply(new CoefficientImpl(num));
                substSet.add(s);
            }
            ComponentSubstParameters compSubstParams = new ComponentSubstParameters(e, substSet, myParams.mfp);
            argument = compSubst.apply(argument, compSubstParams);
        }
        return argument;
    }

    public static class ExchangeSplitParameters implements Parameters {
        public final MulFlexibleParameters mfp;
        public final byte maxExchangeNodes;
        public final char eBond, excBond;
        public final boolean negativeExchange;

        /**
         * @param mfp MulFlexibleParameters (underlying operations require multiplication)
         * @param maxExchangeNodes the maximum exchange group size
         * @param eBond color of e-bonds
         * @param excBond color to use for exchange e-bonds
         * @param negativeExchange true if odd exchanges should be negative
         */
        public ExchangeSplitParameters(MulFlexibleParameters mfp, byte maxExchangeNodes, char eBond, char excBond, boolean negativeExchange) {
            this.mfp = mfp;
            this.maxExchangeNodes = maxExchangeNodes;
            this.eBond = eBond;
            this.excBond = excBond;
            this.negativeExchange = negativeExchange;
        }
    }
}

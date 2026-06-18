/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.*;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.*;
import etomica.graph.operations.AllIsomorphs.AllIsomorphsParameters;
import etomica.graph.operations.IsoFree.IsoFreeParams;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.SplitOneBiconnected.SplitOneParametersBC;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.Property;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.virial.cluster.VirialDiagrams;

import java.util.*;

public class PYGenerator extends IEGenerator {

    /**
     * Returns the graphs (at each order, up to n+2) in the PY approximation to c.
     * Diagrams will include only f-bonds, using the specified color.
     */
    public static List<Set<Graph>> PYGenerate(int n, char fBond){
        
        Set<Graph> c0 = new GraphList();
        Set<Graph> t0 = new GraphList();
        Set<Graph> h0 = new GraphList();
        //create a Map
        ArrayList<Set<Graph>> cList = new ArrayList<Set<Graph>>();
        ArrayList<Set<Graph>> tList = new ArrayList<Set<Graph>>();
        ArrayList<Set<Graph>> hList = new ArrayList<Set<Graph>>();
        
        //fbond
        Graph fGraph = GraphFactory.createGraph((byte)2, BitmapFactory.createBitmap((byte)2,true));
        fGraph.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
        fGraph.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
        fGraph.getEdge((byte)0, (byte)1).setColor(fBond);
                   
        //create a list  
        c0.add(fGraph);
        h0.add(fGraph);
        
        //populate the Map  
        cList.add(c0);
        tList.add(t0);
        hList.add(h0);
        //Clearing the sets since it has been added
        
        IsoFree isoFree = new IsoFree();
        AllIsomorphs allIso = new AllIsomorphs();
        AllIsomorphsParameters aip = new AllIsomorphsParameters(true, new Property() {
            public boolean check(Graph graph) {
                return graph.getNode((byte)0).getType() == Metadata.TYPE_NODE_ROOT && graph.getNode((byte)1).getType() == Metadata.TYPE_NODE_ROOT;
            }
        });
        
        for(int m=1;m<=n;m++){ //Go from C0 to Cn
            Set<Graph> cm = new HashSet<Graph>();
            Set<Graph> tm = new HashSet<Graph>();
            Set<Graph> hm = new HashSet<Graph>();
            
        //  System.out.println(" m = "+m);
            //Calculating tn
            //R
        //  System.out.println(" R ");
            for(int i=0;i<m;i++){
                
                int j=m-i-1;
                
                if(j<0)throw new RuntimeException("oops");
                
                Set<Graph> ci = cList.get(i);
                Set<Graph> hj = hList.get(j);
                //System.out.println("We are multiplying c"+i+" and h"+j+" to get t"+m);
                if(ci.isEmpty()) throw new RuntimeException(i+"cc is empty");
                if(hj.isEmpty()) throw new RuntimeException(j+"hh is empty");
                for(Graph gc:ci){
                    for(Graph gh:hj){
                        //System.out.println("Graph #"+q);q++;
                        Graph newGraph = IEGenerator.fourierMul(gc,gh,(byte)1);
                        newGraph = IEGenerator.relabel(newGraph);
//                        System.out.println("t"+m+" "+newGraph);
                        tm.add(newGraph.copy());
                        hm.add(newGraph.copy());
                    }
                }
                    
            }
            if(tm.isEmpty()) throw new RuntimeException("t set is empty");
            if(hm.isEmpty()) throw new RuntimeException("h set is empty");

            //Calculating c
            for(int i=0;i<m;i++){
            
                int j=m-i-1;        
                
                if(j<0)break;
                
                Set<Graph> ci = cList.get(i);
                Set<Graph> cj = cList.get(j); //kk also contains c
                //System.out.println("We are multiplying c"+i+" and c"+j+" to get c"+m);
                for(Graph gci:ci){
                    for(Graph gcj:cj){
                        
                        //System.out.println("Graph #"+q);q++;
                        Graph temp = IEGenerator.fourierMul(gci,gcj,(byte)1);
                        temp = IEGenerator.relabel(temp);
                        //IEGenerator.GraphDetails(temp);
                        //System.out.println("With the fbond");
                        temp.putEdge((byte)0,(byte)1);
                        temp.getEdge((byte)0,(byte)1).setColor(fBond);
                        //temp = IEGenerator.realMul(temp,fbond).copy();
//                        System.out.println("c"+m+" "+temp);
                        cm.add(temp.copy());
                        hm.add(temp.copy());
                        
                    }
                }
                
                
            }
            //Q
            //System.out.println(" Q ");
            for(int i=0;i<m;i++){
            
                int j=m-i-1;
                    
                if(j<0)break;
                if(j==0) break;//t0 is zero
                
                Set<Graph> ci = cList.get(i);
                Set<Graph> tj = tList.get(j); //kk also contains c
                //System.out.println("We are multiplying c"+i+" and t"+j+" to get c"+m);
                for(Graph gci:ci){
                    for(Graph gtj:tj){
                    
                    //  System.out.println("Graph #"+q);q++;
                        Graph temp = IEGenerator.fourierMul(gci,gtj,(byte)1);
                        temp = IEGenerator.relabel(temp);
                        //System.out.println("With the fbond");
                        temp.putEdge((byte)0,(byte)1);
                        temp.getEdge((byte)0,(byte)1).setColor(fBond);
                        //temp = IEGenerator.realMul(temp,fbond).copy();
                    //  IEGenerator.GraphDetails(temp);
//                        System.out.println("c"+m+" "+temp);
                        cm.add(temp.copy());
                        hm.add(temp.copy());
                    }
                }
                
            }           
            if(cm.isEmpty()) throw new RuntimeException("c set is empty");
            if(hm.isEmpty()) throw new RuntimeException("h set is empty");
            
            //Adding the Graph Set to cmap
            cm = allIso.apply(isoFree.apply(cm, null), aip);
            tm = allIso.apply(isoFree.apply(tm, null), aip);
            hm = allIso.apply(isoFree.apply(hm, null), aip);

            cList.add(cm);
            hList.add(hm);
            tList.add(tm);
        }
        
        
    /*  for(int i=0;i<=n;i++){
            Set<Graph> temp = (Set<Graph>) tmap.get(i);
            System.out.println(" T"+i+" has "+temp.size()+" elements");
        }*/
        
    /*  int p=2;
        Set<Graph> temp = (Set<Graph>) cmap.get(p);
        ClusterViewer.createView("c"+p, temp);*/
    
            
        return cList;
        
    }

    public static Set<Graph> getICPYCorrection(byte n, boolean useE) {
        MetadataImpl.rootPointsSpecial = true;
        final HashMap<Character,Integer> colorOrderMap = VirialDiagrams.initMetaDataComparator();
        if (colorOrderMap != null) {
            colorOrderMap.put('f', 3);
            colorOrderMap.put('e', 4);
        }

        GraphIterator it = new IteratorEF(n, useE);
        IsomorphismFilter iso = new IsomorphismFilter(it);
        MaxIsomorph maxIso = new MaxIsomorph();
        Property root01 = new Property() {
            public boolean check(Graph graph) {
                return graph.getNode((byte)0).getType() == Metadata.TYPE_NODE_ROOT && graph.getNode((byte)1).getType() == Metadata.TYPE_NODE_ROOT;
            }
        };
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), root01);
        
        Set<Graph> pycc = VirialDiagrams.makeGraphList();
        while (iso.hasNext()) {
            Graph g = iso.next();
            Graph giso = maxIso.apply(g, mip);
            pycc.add(giso);
        }
        return pycc;
    }

    /**
     * Returns graphs included in the correction to the PY approximation at the
     * n-th order.  The graphs will be fully-connected with e- and f-bonds.
     */
    public static Set<Graph> getPYCorrection(final byte n) {
        MetadataImpl.rootPointsSpecial = true;
        GraphIterator it = new PropertyFilter(new DefaultIterator(n, (byte)2), new IsBiconnected());

        Set<Graph> cFull = new HashSet<Graph>();
        char fBond = 'f';
        Coefficient fac = new CoefficientImpl(1, (int)SpecialFunctions.factorial(n-2));
        while (it.hasNext()) {
            Graph g = it.next();
            for (Edge e : g.edges()) {
                e.setColor(fBond);
            }
            g.coefficient().multiply(fac);
            cFull.add(g);
        }
//        ClusterViewer.createView("cExact", cFull);
        
        List<Set<Graph>> cList = PYGenerator.PYGenerate(n-2, fBond);
        Set<Graph> cm = cList.get(n-2);

        MulScalar mulScalar = new MulScalar();
        MulScalarParameters msp = new MulScalarParameters(-1, 1);

        Set<Graph> bPYC = new HashSet<Graph>();
        bPYC.addAll(cFull);
        bPYC.addAll(mulScalar.apply(cm, msp));
        IsoFree isoFree = new IsoFree();
        bPYC = isoFree.apply(bPYC, null);
        
        final char nfBond = 'F';
        final char eBond = 'e';
        SplitOneParametersBC splitOneParameters = new SplitOneParametersBC(fBond, eBond, nfBond);
        SplitOneBiconnected splitOneBC = new SplitOneBiconnected();
        Set<Graph> newB = new GraphList();
        for (Graph g : bPYC) {
            Set<Graph> gSet = splitOneBC.apply(g, splitOneParameters);
            for (Graph g2 : gSet) {
                boolean even = true;
                for (Edge e : g2.edges()) {
                    if (e.getColor() == nfBond) {
                        even = !even;
                        e.setColor(fBond);
                    }
                }
                if (!even) {
                    g2 = mulScalar.apply(g2, msp);
                }
                newB.add(g2);
            }
        }
        // all diagrams are fully-connected.  just try to distinguish between them based on e-bonds.
        IsoFreeParams isoFreeParams = new IsoFreeParams() {
            List<Byte> rootMap = new ArrayList<Byte>();
            List<Byte> fieldMap = new ArrayList<Byte>();
            public String getSignature(Graph g) {
                rootMap.clear();
                fieldMap.clear();
                int numE = 0;
                for (byte nodeId = 0; nodeId < n; nodeId++) {
                  byte outDegree = 0;
                  byte nodeCount = n;
                  for (byte i = 0; i < nodeCount; i++) {
                    if (i != nodeId) {
                      outDegree += (g.hasEdge(nodeId, i) && g.getEdge(nodeId, i).getColor() == eBond) ? 1 : 0;
                    }
                  }
                  numE += outDegree;
                  if (nodeId < 2) {
                      rootMap.add(outDegree);
                  }
                  else {
                      if (fieldMap.size() > outDegree) {
                          fieldMap.set(outDegree, (byte)(fieldMap.get(outDegree)+1));
                      }
                      else {
                          for (int i=fieldMap.size(); i<outDegree; i++) {
                              fieldMap.add((byte)0);
                          }
                          fieldMap.add((byte)1);
                      }
                  }
                }
                Collections.sort(rootMap);
                numE /= 2;
                String result = "/E"+numE;
                result += "/R";
                for (Byte degree : rootMap) {
                    result += ":" + degree;
                }
                result += "/F";
                for (Byte degree : fieldMap) {
                    result += ":" + degree;
                }

                return result;
            }
        };
        MaxIsomorph maxIso = new MaxIsomorph();
        Property root01 = new Property() {
            public boolean check(Graph graph) {
                return graph.getNode((byte)0).getType() == Metadata.TYPE_NODE_ROOT && graph.getNode((byte)1).getType() == Metadata.TYPE_NODE_ROOT;
            }
        };
        Set<Graph> newNewB = VirialDiagrams.makeGraphList();
        // take the max isomorph for each graph with the root points still in place
        newNewB.addAll(maxIso.apply(isoFree.apply(newB, isoFreeParams), new MaxIsomorph.MaxIsomorphParameters(new GraphOp.GraphOpNull(), root01)));
        return newNewB;
    }

    public static void main(String[] args){
        MetadataImpl.rootPointsSpecial = true;
        int m = 3;

        if (m < 4) {
            // only try to display diagrams for B6 and below
            List<Set<Graph>> cList = PYGenerator.PYGenerate(m, 'f');
            Set<Graph> cm = VirialDiagrams.makeGraphList();
            cm.addAll(cList.get(m));
            System.out.println("c"+m+"PY");
            for (Graph g : cm) {
                System.out.println(g);
            }
            ClusterViewer.createView("PY", cm);
    
            Set<Graph> correction = VirialDiagrams.makeGraphList();
            correction.addAll(PYGenerator.getPYCorrection((byte)(m+2)));
            System.out.println("c"+(m)+" - c"+(m)+"PY");
            for (Graph g : correction) {
                System.out.println(g);
            }

            ClusterViewer.createView("correctionEF", correction);

            IsoFree isoFree = new IsoFree();
            Set<Graph> correctionB = VirialDiagrams.makeGraphList();
            MulScalar mulScalar = new MulScalar();
            MulScalarParameters msp = new MulScalarParameters(-1, m+2);
            Set<Graph> cb = mulScalar.apply(correction, msp);
            for (Graph g : cb) {
                g.getNode((byte)0).setType(Metadata.TYPE_NODE_FIELD);
                g.getNode((byte)1).setType(Metadata.TYPE_NODE_FIELD);
            }
            MaxIsomorph maxIso = new MaxIsomorph();
            correctionB.addAll(maxIso.apply(isoFree.apply(cb, null), MaxIsomorph.PARAM_ALL));
            ClusterViewer.createView("B correctionEF", correctionB);
        }

        Set<Graph> correctionIC = PYGenerator.getICPYCorrection((byte)(m+2), false);
        System.out.println("c"+(m)+" - c"+(m)+"ICPY (F)");
        for (Graph g : correctionIC) {
            System.out.println(g);
        }
        ClusterViewer.createView("ICcorrectionF", correctionIC);

        correctionIC = PYGenerator.getICPYCorrection((byte)(m+2), true);
        System.out.println("c"+(m)+" - c"+(m)+"ICPY (EF)");
        for (Graph g : correctionIC) {
            System.out.println(g);
        }
        ClusterViewer.createView("ICcorrectionEF", correctionIC);

        IsoFree isoFree = new IsoFree();
        Set<Graph> correctionICB = VirialDiagrams.makeGraphList();
        MulScalar mulScalar = new MulScalar();
        MulScalarParameters msp = new MulScalarParameters(-1, m+2);
        Set<Graph> cb = mulScalar.apply(correctionIC, msp);
        for (Graph g : cb) {
            g.getNode((byte)0).setType(Metadata.TYPE_NODE_FIELD);
            g.getNode((byte)1).setType(Metadata.TYPE_NODE_FIELD);
        }
        MaxIsomorph maxIso = new MaxIsomorph();
        correctionICB.addAll(maxIso.apply(isoFree.apply(cb, null), MaxIsomorph.PARAM_ALL));
        ClusterViewer.createView("B ICcorrectionEF", correctionICB);
    }

    public static class IteratorEF implements GraphIterator {

        private Iterator<Graph> iterator;
        protected final List<Graph> substSet;
        protected final char eBond, fBond, nfBond;
        protected final MulScalar mulScalar;
        protected final MulScalarParameters mspMinus, mspPlus;
        protected final SplitOneParametersBC splitOneParameters;
        protected final SplitOneBiconnected splitOneBC;
        protected final boolean useE;

        public IteratorEF(byte n, boolean useE) {
            this.iterator = new IsomorphismFilter(new PropertyFilter(new PropertyFilter(new DefaultIterator(n, (byte)2), new Property() {
                public boolean check(Graph graph) {
                    return !graph.hasEdge((byte)0, (byte)1);
                }
            }), new IsBiconnected()));
            this.useE = useE;
            eBond = 'e';
            fBond = 'f';
            nfBond = 'F';
            mulScalar = new MulScalar();
            mspMinus = new MulScalarParameters(-1, (int)SpecialFunctions.factorial(n-2));
            mspPlus = new MulScalarParameters(1, (int)SpecialFunctions.factorial(n-2));

            splitOneParameters = new SplitOneParametersBC(fBond, eBond, nfBond);
            splitOneBC = new SplitOneBiconnected();
            substSet = new ArrayList<Graph>();
        }
      
        public boolean hasNext() {
            return substSet.size() > 0 || iterator.hasNext();
        }

        public Graph next() {
            if (substSet.size() > 0) {
                return substSet.remove(substSet.size()-1);
            }
                
            Graph g = iterator.next().copy();
            for (Edge e : g.edges()) {
                e.setColor(fBond);
            }
            g.putEdge((byte)0, (byte)1);
            g.getEdge((byte)0, (byte)1).setColor(eBond);
            if (!useE) {
                return mulScalar.apply(g, mspPlus);
            }

            Set<Graph> gSet = splitOneBC.apply(g, splitOneParameters);
            for (Graph g2 : gSet) {
                boolean even = true;
                for (Edge e : g2.edges()) {
                    if (e.getColor() == nfBond) {
                        even = !even;
                        e.setColor(fBond);
                    }
                }
                if (!even) {
                    g2 = mulScalar.apply(g2, mspMinus);
                }
                else {
                    g2 = mulScalar.apply(g2, mspPlus);
                }
                substSet.add(g2);
            }
            return substSet.remove(substSet.size()-1);
        }
      
        public void remove() {
            throw new RuntimeException("nope");
        }
      }
}

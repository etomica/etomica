package etomica.virial;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Coefficient;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.SplitOneBiconnected;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.SplitOneBiconnected.SplitOneParametersBC;
import etomica.graph.viewer.ClusterViewer;
import etomica.virial.cluster.VirialDiagrams;

public class PYGenerator extends IEGenerator {

    /**
     * Returns the graphs (at each order, up to n+2) in the PY approximation to c.
     * Diagrams will include only f-bonds, using the specified color.
     */
    public static List<Set<Graph>> PYGenerate(int n, char fBond){
        
        Set<Graph> c0 = new GraphList<Graph>();
        Set<Graph> t0 = new GraphList<Graph>();
        Set<Graph> h0 = new GraphList<Graph>();
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
        
        for(int m=1;m<=n;m++){ //Go from C0 to Cn
            HashSet<Graph> cm = new HashSet<Graph>();
            HashSet<Graph> tm = new HashSet<Graph>();
            HashSet<Graph> hm = new HashSet<Graph>();
            
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
                        if(gc.nodeCount()==0) throw new RuntimeException("gc is empty");
                        if(gh.nodeCount()==0) throw new RuntimeException("gh is empty");
                        Graph newGraph = IEGenerator.fourierMul(gc,gh,(byte)1);
                        newGraph = IEGenerator.relabel(newGraph);
//                        System.out.println("t"+m+" "+newGraph);
                        tm.add(newGraph);
                        hm.add(newGraph);
                            
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
                        temp = IEGenerator.relabel(temp).copy();
                        //IEGenerator.GraphDetails(temp);
                        //System.out.println("With the fbond");
                        temp.putEdge((byte)0,(byte)1);
                        temp.getEdge((byte)0,(byte)1).setColor(fBond);
                        //temp = IEGenerator.realMul(temp,fbond).copy();
//                        System.out.println("c"+m+" "+temp);
                        cm.add(temp);
                        hm.add(temp);
                        
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
                        temp = IEGenerator.relabel(temp).copy();
                        //System.out.println("With the fbond");
                        temp.putEdge((byte)0,(byte)1);
                        temp.getEdge((byte)0,(byte)1).setColor(fBond);
                        //temp = IEGenerator.realMul(temp,fbond).copy();
                    //  IEGenerator.GraphDetails(temp);
//                        System.out.println("c"+m+" "+temp);
                        cm.add(temp);
                        hm.add(temp);
                        
                    }
                }
                
            }           
            if(cm.isEmpty()) throw new RuntimeException("c set is empty");
            if(hm.isEmpty()) throw new RuntimeException("h set is empty");
            
            //Adding the Graph Set to cmap
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

    /**
     * Returns graphs included in the correction to the PY approximation at the
     * n-th order.  The graphs will be fully-connect with e- and f-bonds.
     */
    public static Set<Graph> getPYCorrection(int n) {
        VirialDiagrams diagrams = new VirialDiagrams(n, false, false);
        diagrams.setAllPermutations(false);
        diagrams.setDoReeHoover(false);
        diagrams.setDoShortcut(true);
        diagrams.makeVirialDiagrams();
        Set<Graph> bFull = diagrams.getMSMCGraphs(true, false);

        char fBond = diagrams.fBond;
        
        List<Set<Graph>> cList = PYGenerator.PYGenerate(n-2, fBond);
        Set<Graph> cm = cList.get(n-2);
        Coefficient fac = new CoefficientImpl(-1, n);
        for (Graph g : cm) {
            g.getNode((byte)0).setType(Metadata.TYPE_NODE_FIELD);
            g.getNode((byte)1).setType(Metadata.TYPE_NODE_FIELD);
            g.coefficient().multiply(fac);
            System.out.println(g);
        }
        MaxIsomorph maxIso = new MaxIsomorph();
        MaxIsomorphParameters mip = MaxIsomorph.PARAM_ALL;
        cm = maxIso.apply(cm, mip);
        IsoFree isoFree = new IsoFree();
        cm = isoFree.apply(cm, null);
        
        
        MulScalar mulScalar = new MulScalar();
        MulScalarParameters msp = new MulScalarParameters(-1, 1);

        
        Set<Graph> bPYC = new HashSet<Graph>();
        bPYC.addAll(bFull);
        bPYC.addAll(mulScalar.apply(cm, msp));
        bPYC = isoFree.apply(bPYC, null);
        
        char nfBond = 'F';
        char eBond = 'e';
        SplitOneParametersBC splitOneParameters = new SplitOneParametersBC(fBond, eBond, nfBond);
        SplitOneBiconnected splitOneBC = new SplitOneBiconnected();
        Set<Graph> newB = new GraphList<Graph>();
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
        return isoFree.apply(newB, null);
    }

    public static void main(String[] args){

        int m = 3;
        
        List<Set<Graph>> cList = PYGenerator.PYGenerate(m, 'f');
        Set<Graph> cm = cList.get(m);
        Coefficient fac = new CoefficientImpl(-1, m+2);
        for (Graph g : cm) {
            g.getNode((byte)0).setType(Metadata.TYPE_NODE_FIELD);
            g.getNode((byte)1).setType(Metadata.TYPE_NODE_FIELD);
            g.coefficient().multiply(fac);
            System.out.println(g);
        }
        MaxIsomorph maxIso = new MaxIsomorph();
        MaxIsomorphParameters mip = MaxIsomorph.PARAM_ALL;
        cm = maxIso.apply(cm, mip);
        IsoFree isoFree = new IsoFree();
        cm = isoFree.apply(cm, null);
        ClusterViewer.createView("PY", cm);
        
        Set<Graph> correction = PYGenerator.getPYCorrection(m+2);
        ClusterViewer.createView("correctionEF", correction);
    }

}

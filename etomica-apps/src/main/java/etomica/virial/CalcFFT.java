/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;


import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.comparators.*;
import etomica.graph.operations.GraphOp;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.property.FFTDecomposition;
import etomica.graph.property.IsBiconnected;
import etomica.math.SpecialFunctions;
import etomica.math.function.IFunction;
import etomica.math.numerical.SineTransform;
import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * CalcFFT can calculate distribution functions using FFT.
 * Results are cached and reused if the component is needed again.
 * Numerical results for components can also be taken from a file
 *
 * @author Andrew Schultz
 */
public class CalcFFT {

	protected double reducedTemp; // kT/epsilon
	
	protected int N;

	protected SineTransform dst;
    protected double del_r;
	
	protected double[] fr;
	protected double[] fk;
	
	protected double[] ffr;
	protected final HashMap<Object,double[]> rDistributionHash;
    protected final HashMap<Object,double[]> rDistributionErrorHash;
    protected final HashMap<Object,double[]> kDistributionHash;

    public static IFunction makeF(final IPotential2 p2, final double temperature) {
        IFunction f = new IFunction() {
            
            public double f(double r) {
                double r2 = r*r;
                double u = p2.u(r2);
                double x = -u/temperature;
                
                if ( Math.abs(x) < 0.01) {
                    return x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
                }
                return Math.exp(x)-1.0;
            }
        };
        return f;
    }
    
	public CalcFFT(IFunction f, double dr, int power){
		reducedTemp = 1.0;
		N = 1<<power;
		del_r = dr;

	    rDistributionHash = new HashMap<Object,double[]>();
        rDistributionErrorHash = new HashMap<Object,double[]>();
        getfr(f);
        rDistributionHash.put(2, fr);

        kDistributionHash = new HashMap<Object,double[]>();
        dst = new SineTransform();
        fk = dst.forward(fr, del_r);
        kDistributionHash.put(2, fk);

        double[] ffk = new double[N];
        for (int i = 0;i<N; i++) {
            ffk[i] = fk[i]*fk[i];
        }
        kDistributionHash.put(3, ffk);
        
        ffr = dst.reverse(ffk, del_r);
        rDistributionHash.put(3, ffr);
	}

	protected void getfr(IFunction f) {
		
		fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
            double r = n*del_r;
		    fr[n] = f.f(r);
		}
		
	}
	
	public static void Mul(double[] a, double[] b){
		for (int i = 0;i<a.length; i++) {
			a[i] = a[i]*b[i];
		}
	}

	public static void printStrands(List list, int deep) {
	    
		for (Object obj : list) {
		
	      if (obj instanceof Integer) {

	        System.out.print(obj);
	      
	      } 				
	      else {
	       
	        System.out.print("*");
	        System.out.print("(");
	        printStrands((List)obj, deep+1);
	        System.out.print(")");
	      }
	      
	    }
	    
	    if (deep==1) System.out.println(")");
	   	    
	}
	
	public double[][] value(List obj,boolean isK){
        double[] rv = isK ? kDistributionHash.get(obj) : rDistributionHash.get(obj);
        if (rv != null) {
            if (!isK) {
                double[] erv = rDistributionErrorHash.get(obj);
                return new double[][]{rv,erv};
            }
            return new double[][]{rv};
        }
	  	 
	    double[] a = new double[N];
	    double[][] e = new double[0][0];
		 
	    for(int i=0;i<N;i++){
	        a[i]=1;
	    }
	    HashMap<Object,Integer> errMap = new HashMap<Object,Integer>();
        HashMap<Object,Integer> countMap = new HashMap<Object,Integer>();
	    for(Object o:obj){
	        Integer count = countMap.get(o);
	        if (count == null) {
                count = 1;
	        }
	        else {
	            count = count + 1;
	        }
            countMap.put(o, count);
	        
	        double[] b = isK ? kDistributionHash.get(o) : rDistributionHash.get(o);
	        if (b != null) {
	            for (int i=0; i<e.length; i++) {
	                Mul(e[i], b);
		        }
		        double[] oe = rDistributionErrorHash.get(o);
		        if (!isK && count == 1 && oe!=null) {
                    errMap.put(o, e.length);
		            e = Arrays.copyOf(e, e.length + 1);
		            e[e.length-1] = new double[N];
		            System.arraycopy(a, 0, e[e.length-1], 0, N);
		            Mul(e[e.length-1], oe);
		        }
                Mul(a, b);
		    }
		    else if (o instanceof Number) {
		        double[] Ovalue = readDistribution((Integer)o, false);
                if (Ovalue == null) {
//                    System.out.println("don't have distribution function for "+o);
                    return null;
                }
                rDistributionHash.put(o, Ovalue);
                double[] evalue = readDistribution((Integer)o, true);
                if (evalue != null) {
                    rDistributionErrorHash.put(o, evalue);
                }
                if (isK) {
                    Ovalue = dst.forward(Ovalue, del_r);
                    kDistributionHash.put(o, Ovalue);
                }
                // we need to deal with error first because we need a (unmultiplied by o value)
                for (int i=0; i<e.length; i++) {
                    Mul(e[i], Ovalue);
                }
                if (!isK && evalue != null && count == 1) {
                    errMap.put(o, e.length);
                    e = Arrays.copyOf(e, e.length+1);
                    e[e.length-1] = new double[N];
                    System.arraycopy(a, 0, e[e.length-1], 0, N);
                    Mul(e[e.length-1], evalue);
                }
                Mul(a,Ovalue);
		    }
		    else {
		        // if isK, we want to a k-space function, but the elements of this
		        // list need to be multiplied together in r-space.  we'll take the
		        // result and transform to k-space.
	      		double[][] Ovalue = value((List) o,!isK);
	      		    
	      		if (Ovalue == null) return null;
	      		if (isK) {
	      		    // we got back a real-space function, perhaps with errors... ignore them!
	                Ovalue[0] = dst.forward(Ovalue[0], del_r);
	      		    kDistributionHash.put(o, Ovalue[0]);
	      		}
	      		else {
                    Ovalue[0] = dst.reverse(Ovalue[0], del_r);
	      		    rDistributionHash.put(o, Ovalue[0]);
	      		}
	      		Mul(a,Ovalue[0]);
                for (int i=0; i<e.length; i++) {
                    Mul(e[i], Ovalue[0]);
                }
                if (!isK && Ovalue.length>1 && count == 1) {
                    errMap.put(o, e.length);
                    e = Arrays.copyOf(e, e.length+1);
                    e[e.length-1] = new double[N];
                    System.arraycopy(a, 0, e[e.length-1], 0, N);
                    Mul(e[e.length-1], Ovalue[1]);
                }
		    }
		  }

		  if(isK) {
		      return new double[][]{a};
		  }


		  double[] eProduct = new double[N];
          for (int i=0; i<N; i++) {
              eProduct[i] = 0;
          }
		  HashSet seenObj = new HashSet();
		  for (Object o : obj) {
		      if (seenObj.contains(o)) continue;
		      seenObj.add(o);
		      Integer errIndex = errMap.get(o);
		      if (errIndex != null) {
		          int count = countMap.get(o);
		          for (int i=0; i<N; i++) {
		              double ei = e[errIndex][i] * count;
		              eProduct[i] += ei*ei;
		          }
		      }
		  }
          for (int i=0; i<N; i++) {
              eProduct[i] = Math.sqrt(eProduct[i]);
          }
		  return new double[][]{a,eProduct};
	}
	 
    public double[] readDistribution(int graphNum, boolean err) {
//        if (true) return null;
        BufferedReader br = null;
        try {
            File file = new File("d"+graphNum+".dat");
            if (!file.exists()) return null;
            FileReader fileReader = new FileReader(file);
            br = new BufferedReader(fileReader);
            String line = null;
            int i = 0;
            double[] d = new double[N];
            while ((line = br.readLine()) != null) {
                String[] split = line.split(" ");
                if (split.length < 2) throw new RuntimeException("oops bad line "+line);
                if (split.length < 3 && err) return null;
                if (i==1) {
                    double thisdr = Double.parseDouble(split[0]);
                    if (Math.abs(thisdr-del_r) > 1e-14*(thisdr+del_r)) {
                        throw new RuntimeException("dr from file does not match our dr");
                    } 
                }
                d[i] = Double.parseDouble(split[err ? 2 : 1]);
                i++;
            }
            return d;
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        finally {
            if (br != null) {
                try {
                    br.close();
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }
	 
     public static void main(String args[]){
		 
		int n = 6;
		if (args.length>0) {
		     n = Integer.parseInt(args[0]);
		}
		int power = 18;
		double dr = 0.01;
		
        Space space = Space3D.getInstance();
        IPotential2 p2 = new P2LennardJones(1, 1);
		CalcFFT calcFFT = new CalcFFT(makeF(p2, 1.0), dr, power);
//		List myList = new ArrayList();
//		List mySubList = new ArrayList();
//		mySubList.add(31);
//		mySubList.add(2);
//		myList.add(mySubList);
//		double[] dist = calcFFT.value(myList, false)[0];
//		for (int i=0; i<dist.length; i++) {
//		    System.out.println(i*dr+" "+dist[i]);
//		}
//		System.exit(1);
		IsBiconnected isBi = new IsBiconnected();
		GraphIterator iter = new PropertyFilter(new StoredIterator((byte)n), isBi);
		ComparatorChain comp = new ComparatorChain();
		comp.addComparator(new ComparatorNumFieldNodes());
		comp.addComparator(new ComparatorBiConnected());
		comp.addComparator(new ComparatorNumEdges());
		comp.addComparator(new ComparatorNumNodes());
	    GraphList graphList = new GraphList(comp);
		MaxIsomorph maxIso = new MaxIsomorph();
		MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), MaxIsomorph.PROPERTY_ALL);
		
		while (iter.hasNext()) {
		   graphList.add(maxIso.apply(iter.next(), mip));
		}
		long t1 = System.currentTimeMillis();
		FFTDecomposition isFFT = new FFTDecomposition();
		Set<Graph> fftSet = new GraphList(comp);
        Set<Graph> notfftSet = new GraphList(comp);
		//    System.out.println("FFT");
		
        double BnValue = 0;
        double BnError = 0;
        int nd = 0;
		for (Graph g : graphList) {
//		    if (g.edgeCount() < 6) continue;
		   if (isFFT.check(g)) {
		 //       System.out.println("AP: "+isFFT.getArticulationPair());
		         // calcFFT.printStrands(isFFT.getStrands(),1);
		      //  System.out.println(isFFT.getStrands());
		       if (isFFT.getSegments().size() == 0) { // || isFFT.getSegments().get(0) < 400 || isFFT.getSegments().get(0) > 900) {
		           isFFT.getSegments().clear();
		           continue;
		       }
			   double[][] fftvalue = calcFFT.value(isFFT.getStrands(), false);
			   if (fftvalue == null) {
                   notfftSet.add(g);
			       System.out.print("couldn't compute "+g.getStore().toNumberString()+" ");
	               System.out.println(isFFT.getStrands());
	               System.out.println(g);
	               isFFT.getSegments().clear();
			       continue;
			   }
			   nd++;
//               System.out.println(isFFT.getStrands());
			   double bSum = 0;
			   double eSum = 0;
			   for (int i=0; i<fftvalue[0].length; i++) {
			       double r = i*dr;
			       bSum += fftvalue[0][i]*r*r*dr;
			       double ei = fftvalue[1][i]*r*r*dr;
			       eSum += ei*ei;
			   }
			   eSum = 4.0*Math.PI*Math.sqrt(eSum);
			   bSum *= 4.0*Math.PI;
			   double c = g.coefficient().getValue();
			   BnValue += bSum*c;
			   BnError += eSum*Math.abs(c);
		        
				fftSet.add(g);
				
				System.out.println(g.getStore().toNumberString());
//				System.out.println(" => "+BnValue);
				isFFT.getSegments().clear();
		    }
		      
		    //  max=0;
		}
		System.out.println("handled "+nd+" diagrams");
        long t2 = System.currentTimeMillis();
        System.out.println("computation time "+(t2-t1)/1000.0);
		  //  ClusterViewer.createView("not FFT", notFftSet);
		    
	    
//	    ClusterViewer.createView("FFT", fftSet);
//        ClusterViewer.createView("Not FFT", notfftSet);
		  double nfact=SpecialFunctions.factorial(n);
		  
		  BnValue *= -(n-1)/nfact;
		  BnError = (n-1)/nfact*BnError;
		  System.out.println(" B"+n+" = "+BnValue+"   error: "+BnError);  
		    
	}
	
	
	
}

package etomica.virial;

/**
 * This class calculates the sum of all chain clusters using an adaptation of Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterChain extends ClusterSinglyConnected {

    protected final double[] f1, f2, f3;
    protected double value;
    protected double beta;
    
    public ClusterChain(int nPoints, MayerFunction f) {
        super(nPoints, f);
        f1 = fL;//new double[nf];
        f2 = fN;//new double[nf];
        f3 = new double[nf];
        for(int i=0; i<n; i++) {
            f1[1<<i] = 1.0;
            f2[1<<i] = 1.0;
            f3[1<<i] = 1.0;
            for(int j=i+1; j<n; j++) {
                f2[(1<<i)|(1<<j)] = 0.0;//2-vertex graphs don't have 2 or 3 bonds to any vertex
                f3[(1<<i)|(1<<j)] = 0.0;
            }
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterChain c = new ClusterChain(n, f);
        c.setTemperature(1/beta);
        return c;
    }

    /*
     * Computation of sum of purely singly-connected diagrams.
     */
    protected void calcValue() {
        
        super.calcValue();
        
        //no loop needed here because f2 points to same array as fN
        // v=1, no restriction on clusters that don't have the 1-vertex
//        for (int i=2; i<nf; i+=2) { // all even sets; none contain 1
//            //fA[i] = 0;
//            f2[i] = (fL[i] + fN[i]);
//        }

        for (int i=3; i<nf; i+=2) { // every set will contain 1
            //f1 is initialized to fL, so it is all filled in correctly at this point
            //f2 is initialized to fN; we need to remove branch points to correct it
            f3[i] = 0.0;  //initialize sum; when done looping this will be sum of graphs where branch is at v=1

            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            if (jBits==i) continue; //only two vertices, nothing left to do
            //jBits has 1 and next lowest bit in i
            int iii = ii ^ iLow2Bit;//i with 2 lowest bits off
            int jInc = (iii & -iii);//3rd lowest bit, also increment for j
            for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i containing jBits
                int jComp = (i & ~j); //subset of i complementing j
                while ((j|jComp) != i && j<i) { //accelerate sweep through bits to find partition of i
                    int jHighBits = j^jBits;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j==i) break;//got all bits before finding a proper subset
                
                f3[i] += f1[j] * (f2[jComp|1] + f3[jComp|1]);
            }
            f2[i] -= f3[i];//remove from f2 graphs having branch at v=1
        }

        for (int v=1; v<n; v++) {
            int vs1 = 1<<v;
            for (int i=vs1+1; i<nf; i++) {
                
                //skip i if it doesn't contain v...
                if ((i & vs1) == 0) continue;//
                //...or if it has only one bit
                int iLowBit = (i&-i);//lowest bit in i
                if (iLowBit == i) continue;//low bit is only bit
                
                //start f1(v) as sum of all graphs with no branch at any vertex less than v
                //then subtract f2(v) and f3(v) from it after loop
                f1[i] += f2[i];//f1(v) = f1(v-1) + f2(v-1)
                
                //we'll sum over partitions to get f2(v) and f3(v)
                f2[i] = 0.0;
                f3[i] = 0.0;

                int jBits;
                int ii = i ^ iLowBit;
                int iLow2Bit = (ii & -ii);
                if (iLowBit != vs1 && iLow2Bit != vs1) {
                    //v is not in the lowest 2 bits
                    // jBits is the lowest bit and v
                    jBits = iLowBit | vs1;

                    // we can only increment by the 2nd lowest
                    int jInc = iLow2Bit;

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                        if ((j & jBits) != jBits) {
                            //ensure jBits are in j
                            j |= vs1;
                            if (j==i) break;
                        }
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j|jComp) != i && j<i) {
                            int jHighBits = j^jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow; // this might knock out the v bit
                            j |= vs1;
                            jComp = (i & ~j);
                        }
                        if (j==i) break;
                        f2[i] += f1[j] * f1[jComp|vs1];
                        f3[i] += f1[j] * (f2[jComp|vs1] + f3[jComp|vs1]);
                    }
                }
                else {
                    //lowest 2 bits contain v
                    // jBits is the lowest 2 bits
                    // we can start at jBits and increment by the 3rd lowest bit
                    jBits = iLowBit | iLow2Bit;
                    if (jBits == i) continue; // no bits left for jComp
                    
                    int iii = ii ^ iLow2Bit;
                    int jInc = (iii & -iii);

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                        // start=jBits and jInc ensure that every set includes jBits
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j|jComp) != i && j<i) {
                            int jHighBits = j^jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow;
                            jComp = (i & ~j);
                        }
                        if (j==i) break;
                        f2[i] += f1[j] * f1[jComp|vs1];
                        f3[i] += f1[j] * (f2[jComp|vs1] + f3[jComp|vs1]);
                    }
                }

                f1[i] -= (f2[i] + f3[i]);//remove from f1 graphs where v is not a leaf
            }
        }

        value = f1[nf-1] + f2[nf-1]; 

    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                f1[(1<<i)|(1<<j)] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
}

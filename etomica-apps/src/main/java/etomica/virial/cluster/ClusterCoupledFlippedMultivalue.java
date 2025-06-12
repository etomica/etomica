/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.histogram.HistogramExpanding;
import etomica.math.DoubleRange;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;

public class ClusterCoupledFlippedMultivalue implements ClusterAbstractMultivalue {

    protected final ClusterAbstractMultivalue wrappedCluster;
    protected final ClusterAbstractMultivalue wrappedClusterBD;
    protected final Space space;
    protected final boolean[] flippedAtoms;
    protected final double minFlipDistance;
    protected final int nDer;
    protected  final double tol;
    protected final boolean countflips = true;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value[], lastValue[];
    protected IMoleculePositionDefinition positionDefinition;
    protected long flipcount = 0;
    protected long BDcount = 0;
    protected long totcount = 0;
    protected boolean pushme=false;
    private Vector childAtomVector;
    public HistogramExpanding histe;
    protected boolean[][] printed = new boolean[100][2];
    public boolean doPrint = false;
    public IRandom random;
    public double BDAccFrac = 1;
    public double FlipAccFrac;
    public boolean valueBD;
    public boolean lastValueBD;
    
    /**
     * cluster must have caching disabled
     * configurations will be flipped when the minimum distance between any two molecules
     * exceeds minFlipDistance.  set minFlipDistance to 0 to always flip.
     */
    public ClusterCoupledFlippedMultivalue(ClusterAbstractMultivalue cluster, ClusterAbstractMultivalue clusterBD,Space space, double minFlipDistance, int nDer, double tol) {
        this.space = space;
        wrappedCluster = cluster;
        wrappedClusterBD = clusterBD;
        childAtomVector = space.makeVector();
        flippedAtoms = new boolean[cluster.pointCount()];
        positionDefinition = new MoleculePositionGeometricCenter(space);
        this.minFlipDistance = minFlipDistance;
        value = new double[nDer+1];
        lastValue = new double[nDer+1];
        this.nDer = nDer;
        this.tol = tol;
        histe = new HistogramExpanding(1.0,new DoubleRange(-40,0));
        FlipAccFrac= 1/Math.sqrt(Math.pow(2, wrappedCluster.pointCount()));
    }

    public ClusterAbstract makeCopy() {
        if (wrappedClusterBD==null){
            return new ClusterCoupledFlippedMultivalue((ClusterAbstractMultivalue)wrappedCluster.makeCopy(),null, space, minFlipDistance, nDer, tol);
        }
        else {
            return new ClusterCoupledFlippedMultivalue((ClusterAbstractMultivalue) wrappedCluster.makeCopy(), (ClusterAbstractMultivalue) wrappedClusterBD.makeCopy(), space, minFlipDistance, nDer, tol);
        }
    }

    public int pointCount() {
        return wrappedCluster.pointCount();
    }

    public ClusterAbstract getSubCluster() {
        return wrappedCluster;
    }

    public void setBDAccFrac(double p, IRandom rng ){ BDAccFrac = p; random = rng;}

    public boolean valueIsBD() {
        return valueBD;
    }

    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        long thisCPairID = cPairs.getID();
        boolean debugme = false;
//        if (debugme)System.out.println("1 "+thisCPairID+" "+cPairID+" "+lastCPairID+" "+value[0]+" "+lastValue[0]+" "/*+f[0].getClass()*/);
        if (thisCPairID == cPairID) {
//            if (debugme)System.out.println("2 "+"clusterSum "+cPairID+" returning RECENT "+value[0]);
            return value[0];
        }
        else if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            for(int m=0;m<value.length;m++){
                value[m] = lastValue[m];
            }
//            if (debugme)System.out.println("3 "+"clusterSum "+cPairID+" returning PREVIOUS RECENT "+lastValue[0]);
            valueBD = lastValueBD;
            return value[0];
        }

        // a new cluster
        lastCPairID = cPairID;
        for(int m=0; m<value.length; m++){
            lastValue[m] = value[m];
        }
        lastValueBD = valueBD;
        cPairID = thisCPairID;

        final int pointCount = wrappedCluster.pointCount();

        boolean flipit = false;
        totcount++;
        double minR2 = minFlipDistance*minFlipDistance;
        double maxR2 = 0;
        for (int i=0; i<pointCount; i++) {
            for (int j=i+1; j<pointCount; j++) {
                if (box.getCPairSet().getr2(i,j)>maxR2) maxR2 = box.getCPairSet().getr2(i,j);
                if (box.getCPairSet().getr2(i,j) > minR2) {
                    if ( box.getCPairSet().getr2(i,j) > 400) debugme=true;
                    flipit=true;
                }
            }
        }
//        if (debugme)System.out.println("FLIPIT " + flipit);
        if(flipit&&countflips) flipcount++;
        if(pushme&&maxR2<minR2){
            value[0] = 1e-20;
            return value[0];
        }

        ClusterAbstractMultivalue currentcluster = wrappedCluster;

        int n = wrappedCluster.pointCount();
        double bfac = (1.0-n)/ SpecialFunctions.factorial(n);
        int printIdx = (int)(Math.sqrt(maxR2)/10);

        valueBD=false;

        double[] varr = new double[(1<<n)+1];
        varr[1] = Double.NaN;
        boolean doBD = random == null|| random.nextDouble()<BDAccFrac;
        boolean doFLip = random == null || random.nextDouble()<FlipAccFrac;

        if ( flipit && !doFLip){
            value[0]=0;
            return value[0];
        }

        for (int pass=0; pass<2; pass++){
            if(pass==1){
                if (wrappedClusterBD == null) {
                    value[0] = 0;
                    return value[0];
                }
                currentcluster = wrappedClusterBD;
            }

            double doubleval = value[0];
            double[] v = currentcluster.getAllLastValues(box);
            if(pass==0)varr[0] = v[0]/bfac;

            if (doPrint && pass==1 && !printed[printIdx][1] ) {
                System.out.println("BD "+v[0]/bfac);
            }
            for (int m = 0; m < value.length; m++) {
                value[m] = v[m];
            }

            if (!flipit) {
                if (Math.abs(value[0]/bfac)<tol && pass==0 && value[0]!=0){
                    if(!doBD){
                        value[0]=0;
                        return value[0];
                    }
                    BDcount++;
                    valueBD=true;
                    if (doPrint && !printed[printIdx][1] ) {
                        varr[1]=Double.NaN;
                        printConfig(box,varr);
                    }
                    continue;
                }
                if( doPrint && !printed[printIdx][pass]){
                    if(pass==0){
                        varr[1]=Double.NaN;
                        printConfig(box,varr);
                    }
                    printed[printIdx][pass] = true;
                }
                return v[0];
            }
            if (debugme) System.out.print(String.format("%6.3f %10.4e ", Math.sqrt(maxR2), v[0]));

            for (int i=0; i<pointCount; i++) {
                flippedAtoms[i] = false;
            }

            IMoleculeList atomList = box.getMoleculeList();
            // loop through the atoms, toggling each one until we toggle one "on"
            // this should generate each combination of flipped/unflipped for all
            // the molecules

            int flipn = 1;
            while (true) {
                boolean didFlipTrue = false;
                for (int i = 0; !didFlipTrue && i < pointCount; i++) {
                    flippedAtoms[i] = !flippedAtoms[i];
                    didFlipTrue = flippedAtoms[i];
                    flip(atomList.get(i));
                }
                if (!didFlipTrue) {
                    // if we flipped every atom from true to false, we must be done
                    break;
                }
                v = currentcluster.getAllLastValues(box);
                if(pass==0){
                    varr[flipn]=v[0]/bfac;
                    flipn++;
                }
                if (doPrint && pass==1 && !printed[printIdx][1]) {
                    System.out.println("BD "+ v[0]/bfac);
                }
                if (debugme) System.out.print(String.format("%10.4e ", v[0]));
                for (int m = 0; m < value.length; m++) {
                    value[m] += v[m];
                }
                if (Double.isNaN(value[0])) {
                    throw new RuntimeException("oops");
                }
            }
            for (int m = 0; m < value.length; m++) {
                value[m] /= (1<<pointCount);
            }
            if(pass==0)varr[varr.length-1] = value[0]/bfac;

            if (doPrint && pass==1 && !printed[printIdx][1]) {
                System.out.println("avg "+value[0]/bfac + "\n");
                printed[printIdx][1] = true;
            }
            if(pass==0&&value[0]!=0)histe.addValue(Math.log(Math.abs(value[0]/bfac)));
            if (debugme) System.out.print(String.format("%10.4e\n", value[0]));

            if(false&&pushme&&pass==1)System.out.println(Math.sqrt(maxR2)+" "+ doubleval+" "+value[0]);

            if ( Math.abs(value[0]/bfac)>tol || value[0]==0) {
                if(doPrint && pass==0 && !printed[printIdx][0]){
                    printConfig(box,varr);
                    printed[printIdx][0]=true;
                }
                break;
            }
            if(pass==0 ){
                if(!doBD){
                    value[0]=0;
                    return value[0];
                }
                BDcount+=1<<pointCount;
                valueBD=true;
                if(doPrint && !printed[printIdx][1]) {
                    printConfig(box, varr);
                }
            }
        }
        double a = (flipit ? FlipAccFrac : 1) * (currentcluster==wrappedClusterBD ? BDAccFrac : 1);

        for (int m = 0; m < value.length; m++) {
            value[m] /= a;
        }
        return value[0];
    }

    private void printConfig(BoxCluster box, double[] varr) {
        System.out.println();
        int nmols = box.getMoleculeList().size();
        for(int i=0; i<nmols;i++){
            IMolecule mol = box.getMoleculeList().get(i);
            int natoms = mol.getChildList().size();
            for(int j=0; j<natoms;j++){
                IAtom atom = mol.getChildList().get(j);
                Vector p = atom.getPosition();
                System.out.println(atom.getType().getElement().getSymbol() + " " + p.getX(0) + " " + p.getX(1) +" "+ p.getX(2));
            }

        }
        for(int i =0; i<varr.length; i++) {
            if(Double.isNaN(varr[i])) break;
            if(i==varr.length-1)System.out.print("avg ");
            System.out.println(varr[i]);
        }
    }

    public int getNumValues() {
        return value.length;
    }

    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }

    protected void flip(IMolecule flippedMolecule) {
        Vector COM = positionDefinition.position(flippedMolecule);
		IAtomList childAtoms = flippedMolecule.getChildList();
		for (int i = 0; i < childAtoms.size(); i++) {
		    childAtomVector.Ea1Tv1(2,COM);
			childAtomVector.ME(childAtoms.get(i).getPosition());
			childAtoms.get(i).getPosition().E(childAtomVector);
		}
    }

    public void setTemperature(double temperature) {
        wrappedCluster.setTemperature(temperature);
        if(wrappedClusterBD!=null)wrappedClusterBD.setTemperature(temperature);
    }

    public long getflipcount(){
        return flipcount;
    }

    public long gettotcount(){
        return totcount;
    }

    public double getflipfrac(){
        return ((double)flipcount)/totcount;
    }

    public long getBDcount() {return BDcount;}

    public long getBDtotcount(){return (totcount-flipcount)+(flipcount*(1<<pointCount()));}

    public double getBDfrac() {return ((double)BDcount)/getBDtotcount();}

}

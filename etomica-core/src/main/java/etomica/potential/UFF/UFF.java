package etomica.potential.UFF;

import etomica.UFF.PDBBuilder;
import etomica.potential.P2Electrostatic;
import etomica.potential.P2Harmonic;
import etomica.potential.P2LennardJones;
import etomica.species.ISpecies;
import etomica.units.Degree;
import etomica.units.Unit;

import java.util.*;

public class UFF {
    static ISpecies species;
    public static double getbondUFF(double ri, double rj, double n, double lambda, double zi, double zj, double chiI, double chiJ){
        double sqrt = Math.pow((Math.sqrt(chiI)-Math.sqrt(chiJ)),2);
        double ren = ri*rj*sqrt/(chiI*ri+chiJ*rj);
        double rbo = - lambda * (ri+rj) * Math.log(n);

        double rjk = ri + rj + rbo + ren;
        return rjk;
    }

    public static P2Harmonic bondUFF(double ri, double rj, double n, double lambda, double zi, double zj, double chiI, double chiJ){
        double rjk = getbondUFF( ri, rj, n,  lambda, zi,   zj,  chiI,  chiJ);
        double kjk = 664.32*zi*zj/Math.pow(rjk,3);
        return new P2Harmonic(kjk, rjk);
    }

    public static P3BondAngleUFF angleUFF(double ri, double rj, double rk, double n, double lambda, double zi, double zj, double zk, double chiI, double chiJ, double chiK, double theta0){
        double rij, rjk, rik, beta, kijk, eTheta, c0, c1, c2, costheta0, cos2theta0;
        costheta0 = Math.cos(theta0);
        cos2theta0 = 2*costheta0*costheta0-1;
        c2=1/Math.pow(4*Math.sin(theta0),2);
        c1=-4*c2*costheta0;
        c0=c2*((2*costheta0*costheta0)+1);
        rij = getbondUFF(ri, rj, n, lambda, zi, zj, chiI, chiJ);
        rjk = getbondUFF(rj, rk, n, lambda, zj, zk, chiJ, chiK);
        rik = Math.sqrt(rij*rij + rjk*rjk -2*rij*rjk*costheta0); //from openBabel
        beta= 664.12*(zi*zk)/Math.pow(rik,5);
        kijk= beta*((3*rij*rjk*(1-costheta0*costheta0))+(rik*rik*costheta0));
        return new P3BondAngleUFF(c0, c1, c2, kijk);
    }


    public static P4BondTorsionUFF torsionUFF(double V, int type){
        double phi = 0;
        Unit degree = Degree.UNIT;

        if(type == 2){
            phi =0;
            double thetha0Rad = degree.toSim(phi);
            return new P4BondTorsionUFF(V,6, thetha0Rad);
        } else if (type ==3) {
            phi =0;
            double thetha0Rad = degree.toSim(phi);
            double bondOrder = 0.7;
            double Vsp2 = 5 * V * ( 1 + 4.18 * Math.log( bondOrder));
            return new P4BondTorsionUFF(Vsp2, 6, thetha0Rad);
        } else {
            phi =0;
            double thetha0Rad = degree.toSim(phi);
            double bondOrder = 0.7;
            double Vsp2 = 5 * V* ( 1 + 4.18 * Math.log( bondOrder));
            return new P4BondTorsionUFF(Vsp2, 6, thetha0Rad);
        }
    }

/*
    public static void torsionTypeIdentifier(int type, double V){
        if(type==2){
            torsionValuesp3(V);
        } else if (type ==3) {
            torsionValueOnesp2(V);
        } else{
            torsionValueTwosp2(V);
        }
    }

*//*
    public static double torsiontype (int i, String atomName){
        double V = 0, Vi = 0, Vj= 0;
        if(i == 1){
            Vi=switchCaseTorsion(atomName);
            System.out.println( Vi + atomName);
        } else if (i == 2) {
            Vj=switchCaseTorsion(atomName);
            System.out.println( Vj + atomName);
        }
        if( Vi > 0 && Vj > 0){
            V = Vi * Vj;
            System.out.println(" V" + V);
        }
        return V;
    }

*/

    public static double switchCaseTorsion (String atomName){
        double V;
        //System.out.println(atomName + "AtomName");
        if(atomName.equals("C")){
            V = 2.119;
            //System.out.println("C " + 2.119 );
        } else if (atomName.equals("O")) {
            V = 0.018;
            //System.out.println("O " + V);
        } else if (atomName.equals("N")){
            V = 0.450;
            //System.out.println("N " + V );
        }else if (atomName.equals("P")){
            V = 1.225;
            //System.out.println("P " + V );
        }else if (atomName.equals("Si")){
            V =  2.400;
            //System.out.println("Si " + V );
        }else {
            V = 2.400;
            //System.out.println("none " + V );
        }
        /*
        switch (atomName){
            case "C":
                System.out.println("C" + 2.119);
                V = 2.119;
            case "N":
                System.out.println("N" + 0.450);
                V = 2.119;
            case "O":
                System.out.println("O" + 0.018);
                V = 2.119;
            case "P":
                System.out.println("P" + 1.225);
                V = 2.119;
            case "Si":
                System.out.println("Si" + 2.400);
                V = 2.119;
            default:
                System.out.println("None" + 2.400);
                V = 0;
        }*/
        return V;
    }
/*species= PDBBuilder.buildSpecies(confName);
        IMolecule molecule = species.makeMolecule();
        ArrayList<ArrayList<Integer>>connectivity=  PDBBuilder.getConnectivity(confName);
        Map<Integer, String> atomMap = PDBBuilder.getAtomMap(connectivity);
        PDBBuilder.bonding(connectivity, atomMap);
        /*for(int i=0; i<torsionArray.length; i++){
            if(i ==2 || i==3 ){

            }
        }*/
/*

    public P4BondInversionUFF inversionUFF (){
        switch (str){
            case invIdenMeth :
                return new P4BondInversionUFF();
            case carbonaromaticChecker:
                return new P4BondInversionUFF(1, -1, 0, 6 );
            case invIdenCarboxyl:
                return new P4BondInversionUFF(1, -1, 0, 50);
            default:
                return new P4BondInversionUFF(0,0,0,0);
        }
    }




 */


    public static void electrostaticUFF(){
        ArrayList<Integer> bondList = new ArrayList<>();
        String confName = "F:/ethene";
        bondList = PDBBuilder.getBonding(confName);
        ArrayList<ArrayList<Integer>> connectivity = PDBBuilder.getConnectivity(confName);
        ArrayList<ArrayList<Integer>> connectivityModified = PDBBuilder.getconnectivityModified(connectivity);
        Map<Integer,String> atomMap = PDBBuilder.getAtomMap(connectivity);
        HashMap<Integer, String> atomMapModified = PDBBuilder.getatomMapModified(atomMap);
        List<int[]> listOfTorsions = PDBBuilder.getTorsionList(connectivity);
        System.out.println(bondList);
        System.out.println(connectivity);
        System.out.println(connectivityModified);
        System.out.println(atomMap);
        System.out.println(atomMapModified);

    }

    public P2Electrostatic electroUFF(double q1, double q2){
        P2Electrostatic p2electro = new P2Electrostatic(q1, q2);
        p2electro.setCharge1(332.0637*q1);
        p2electro.setCharge2(q2);
        return new P2Electrostatic(q1, q2);
    }
    public LJUFF vanderWaals(double xa, double xb, double da, double db, double sciA, double sciB){
        double xab = (xa + xb)/2;
        double dab = da * db;
        double sciAB = (sciA + sciB)/2;
        return new LJUFF(xab, dab, sciAB);
    }

    public P2LennardJones vdw(double xa, double xb, double da, double db){
        double xab = (xa + xb)/2;
        double dab = da * db;

        return new P2LennardJones( xab/Math.pow(2, 1.0/6.0),dab);
    }
}

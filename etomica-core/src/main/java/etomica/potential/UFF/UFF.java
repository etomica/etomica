package etomica.potential.UFF;

import etomica.potential.P2Electrostatic;
import etomica.species.ISpecies;
import etomica.units.*;
import g3dsys.images.Box;

public class UFF {
    static ISpecies species;
    static double[] bondConstantArray;
    public static double getbondUFF(double ri, double rj, double chiI, double chiJ, double bo){
        double part = (Math.sqrt(chiI)-Math.sqrt(chiJ));
        double sqSqrt = part * part;
        double ren = (ri*rj*sqSqrt)/((chiI*ri)+(chiJ*rj));
        double lambda = 0.1332;
        double rbo = - lambda * (ri+rj) * Math.log(bo);
        //double rjk = ri + rj +rbo - ren;
        //System.out.println(bo+ " " + rjk);
        return ri + rj +rbo - ren;
    }

    public static double[] bondUFF(double ri, double rj, double zi, double zj, double chiI, double chiJ, double bondOrder){
        double rjk = getbondUFF( ri, rj,  chiI,  chiJ, bondOrder);
        double kjk = (0.5*664.12*zi*zj/((Math.pow(rjk,3))));//UFF error. Difference of factor of 2.
        Unit[] newOnw = {new PrefixedUnit(Prefix.KILO,Calorie.UNIT), Mole.UNIT, Angstrom.UNIT};
        double[] newExpo ={1.0, -1.0, -2.0};
        CompoundUnit molA2 = new CompoundUnit(newOnw,newExpo);
        double newkijk = molA2.toSim(kjk);
       // double newkijkkJ = kjk * 4.182;
        return new double[]{newkijk, rjk};
    }
    public static double[] BondConstantArray (double kijk, double rjk){
        Unit[] newOnw = {new PrefixedUnit(Prefix.KILO,Joule.UNIT), Mole.UNIT, Angstrom.UNIT};
        double[] newExpo ={1.0, -1.0, -2.0};
        CompoundUnit molA2 = new CompoundUnit(newOnw,newExpo);
        double newkijk = molA2.fromSim(kijk);
        return new double[]{newkijk, rjk};
    }


    public static double valueIdentifier(double atomNumOne, double atomNumTwo){
        double bo =0;
        if( atomNumOne == 1 && atomNumTwo == 1){
            bo = 1;
        } else if (atomNumOne >1 && atomNumTwo ==1 || atomNumOne ==1 && atomNumTwo >1 ) {
            bo = 1;
        } else if (atomNumOne ==2 && atomNumTwo ==2) {
            bo = 2;
        } else if (atomNumOne == 1.5 && atomNumTwo == 1.5) {
            bo = 1.5;
        }  else {
            //bo value for single bond is mentioned as 1
            bo = 1;
        }
      //  System.out.println(bo + " bo");
        return bo;
    }
    public static double[] angleUFF(double ri, double rj, double rk, double zi, double zj, double zk, double chiI, double chiJ, double chiK, double theta0, double valueOne, double valueTwo, double valueThree, int num ){
        double rij, rjk, rik, beta, kijk, c0, c1, c2, costheta0, sintheta0, kb;
        double thetanormal = theta0;
        theta0 = Degree.UNIT.toSim(theta0);
        costheta0 = Math.cos(theta0);
        sintheta0 = Math.sin(theta0);
        c2=1/(4*sintheta0*sintheta0);
        c1=-4*c2*costheta0;
        c0=c2*(2*costheta0*costheta0+1);
        double bo1 = valueIdentifier(valueOne, valueTwo);
        double bo2 = valueIdentifier(valueTwo, valueThree);
        rij = getbondUFF(ri, rj, chiI, chiJ, bo1);
        rjk = getbondUFF(rj, rk, chiJ, chiK, bo2);
        //System.out.println(rij + " " + rjk);
        rik = Math.sqrt(rij*rij + rjk*rjk -2*rij*rjk*costheta0); //from openBabel
        beta = 664.12 / (rij * rjk);
        kijk = beta * (zi*zk/Math.pow(rik,5)) * ((3*rij*rjk*(1-costheta0*costheta0))-(rik*rik*costheta0))* rij*rjk;
        if( num != 0){
            kijk = kijk / 9;
        }
        Unit[] newOnw = {new PrefixedUnit(Prefix.KILO,Calorie.UNIT), Mole.UNIT, Radian.UNIT};
        double[] newExpo ={1.0, -1.0, -2.0};
        CompoundUnit molrad2 = new CompoundUnit(newOnw,newExpo);
        double newkijkkJ = kijk * 4.182;
        double newkijk = molrad2.toSim(kijk);
       // System.out.println(newkijkkJ + " kij");
        return new double[]{(float) c0, (float) c1, (float) c2, (float) newkijk, num};
    }

    public static double[] torsionUFF(double V, int type,double rbo){
        double phi = 0;
        double n;
        Unit degree = Degree.UNIT;
        if(type == 2){
            phi = 60;
            n = 3;
            double thetha0Rad = degree.toSim(phi);
            return new double[]{V, n, thetha0Rad};
        } else if (type == 3) {
            phi = 0;
            n = 6;
            double thetha0Rad = degree.toSim(phi);
            return new double[]{V, n, thetha0Rad};
        } else {
            phi =180;
            n = 2;
            rbo = 1.5;
            double thetha0Rad = degree.toSim(phi);
            double Vsp3 = 5 * V * ( 1 + 4.18 * Math.log(rbo));
            return new double[]{Vsp3, n, thetha0Rad};
        }
    }

    public static double switchCaseTorsion (String atomName){
        double V;
        //System.out.println(atomName + "AtomName");
        if(atomName.equals("C_3")){
            V = 2.119;
            //System.out.println("C " + 2.119 );
        } else if (atomName.equals("O_3")) {
            V = 0.018;
            //System.out.println("O " + V);
        } else if (atomName.equals("N_3")){
            V = 0.450;
            //System.out.println("N " + V );
        }else if (atomName.equals("P_3")){
            V = 1.225;
            //System.out.println("P " + V );
        }else if (atomName.equals("Si")){
            V =  2.400;
            //System.out.println("Si " + V );
        }else {
            V = 2.400;
            //System.out.println("none " + V );
        }
        return V;
    }


    public P2Electrostatic electroUFF(double q1, double q2){
        P2Electrostatic p2electro = new P2Electrostatic(q1, q2);
        p2electro.setCharge1(q1);
        p2electro.setCharge2(q2);
        return new P2Electrostatic(q1, q2);
    }

    public LJUFF vdw(double xa, double xb, double da, double db){
        double xab = (xa + xb)/2;
       // System.out.println(xa + " " + xb + " " + da + " " + db);
        //double xab = Math.sqrt(xa * xb);
        double dab = Math.sqrt(da * db);
        return new LJUFF(xab, dab);
    }


  /*  public LJTraPPE vdwTraPPE(double xa, double xb, double da, double db){
        double xab = Math.sqrt(xa * xb);
        double dab = Math.sqrt(da * db);
        return new LJTraPPE(xab, dab);
    }
*/

}

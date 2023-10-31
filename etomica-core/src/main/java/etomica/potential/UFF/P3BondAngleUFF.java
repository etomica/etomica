package etomica.potential.UFF;

import etomica.potential.IPotentialBondAngle;
import etomica.units.*;
import etomica.util.MathTiming;

public class P3BondAngleUFF  implements IPotentialBondAngle{
    public double theta0;
    private int num, caseNum;
    protected double c0, c1, c2, kijk;
    public P3BondAngleUFF( double c0, double c1, double c2, double kijk, double theta0,  int num, int caseNum) {
        this.c0 = c0;
        this.c1 = c1;
        this.c2 = c2;
        this.kijk = kijk;
        this.num = num;
        this.caseNum = caseNum;
        this.theta0 = theta0;
        //System.out.println(num + " " + caseNum + " " + kijk);
    }

    public double u(double costheta) {
        double u = 0;

//        System.out.println(c0 + " " +c1+ " " + c2 + " " + num + " "+ theta0  + " values" + costheta + " costheta") ;
        if(caseNum == 0){
            double cos2theta = 2 * costheta * costheta - 1;
            u = kijk * ( c0 + c1 * costheta + c2 * cos2theta);
        } else {
            double theta= Math.acos(costheta);
            u = kijk * (1 - Math.cos(num*theta)) + Math.exp(-20.0*(theta - theta0 + 0.25));
           // System.out.println(kijk + " " + (1 - Math.cos(num*theta)) + " "+ Math.exp(-20.0*(theta - theta0 + 0.25)));
        }
        //System.out.println(u + " Inside" );

        return u;
    }

    public void udu(double costheta, double[] u, double[] du) {
        Unit kjmol = new UnitRatio(new PrefixedUnit(Prefix.KILO,Joule.UNIT), Mole.UNIT);
        double theta = Math.acos(costheta);
        double t = theta * 180 / Math.PI;
        double theta0Act = theta0 * Math.PI / 180;
       // theta0 = theta0 * Math.PI / 180;
        double val = theta - theta0Act + 0.25;
        //System.exit(1);
        //System.out.println(theta + " theta here" + theta0 + " theta0 here");
        if(caseNum == 0){
            u[0] =  kijk*  ( c0 + c1 * costheta + c2 * Math.cos(2 * theta));
            du[0] =  -kijk*(c1 + 4*costheta*c2);
           // System.out.println(u[0] + " case 0 "  + kjmol.fromSim(u[0]));
       } else {
            double funct = Math.exp(-20.0*(theta - theta0Act + 0.25));
            //System.out.println(kjmol.fromSim(kijk) + " " + num + " " + theta + " " + theta0Act);
            //double part1 = kjmol.fromSim(kijk * (1 - Math.cos(num * theta)));// double part2 = (Math.exp(-20.0*(theta - theta0Act + 0.25)));
//            double expo =  -20.0*(theta - theta0Act + 0.25);
            u[0] = kijk * (1 - Math.cos(num * theta)) + kjmol.toSim(funct);
            //System.out.println(part1 + " " + part2 + " "+ expo +" "+ kjmol.fromSim(u[0]));
            du[0] = kijk * Math.sin(num * theta) - 20.0 * funct;
            //System.out.println(  kjmol.fromSim(u[0]));
        }
       // System.out.println( " Energy Angle : " + u[0] );
        /*double cosntheta = Math.cos(num*theta);
        u[0] = (kijk / (num *  num))*(1-cosntheta);
        du[0] = -(kijk/num)*(Math.sin(num*theta));*/
    }
}

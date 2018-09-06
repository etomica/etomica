package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;

public class FunctionPhiargpt5parabolicpsinesum  implements FunctionDifferentiable {


    private final double L;
    private final double sa1;
    private final double sb1;
    private final double sc1;
    private final double sa2;
    private final double sb2;
    private final double sc2;
    private final double sa3;
    private final double sb3;
    private final double sc3;
    private final double sa4;
    private final double sb4;
    private final double sc4;
    private final double sa5;
    private final double sb5;
    private final double sc5;
    private final double sa6;
    private final double sb6;
    private final double sc6;
    private final double sa7;
    private final double sb7;
    private final double sc7;
    private final double sa8;
    private final double sb8;
    private final double sc8;

    public FunctionPhiargpt5parabolicpsinesum(double L, double sa1, double sb1, double sc1,double sa2, double sb2, double sc2,double sa3, double sb3, double sc3,double sa4, double sb4, double sc4, double sa5, double sb5, double sc5,double sa6, double sb6, double sc6,double sa7, double sb7, double sc7,double sa8, double sb8, double sc8) {
        this.L = L;
        this.sa1 = sa1;
        this.sb1 = sb1;
        this.sc1 = sc1;
        this.sa2 = sa2;
        this.sb2 = sb2;
        this.sc2 = sc2;
        this.sa3 = sa3;
        this.sb3 = sb3;
        this.sc3 = sc3;
        this.sa4 = sa4;
        this.sb4 = sb4;
        this.sc4 = sc4;
        this.sa5 = sa5;
        this.sb5 = sb5;
        this.sc5 = sc5;
        this.sa6 = sa6;
        this.sb6 = sb6;
        this.sc6 = sc6;
        this.sa7 = sa7;
        this.sb7 = sb7;
        this.sc7 = sc7;
        this.sa8 = sa8;
        this.sb8 = sb8;
        this.sc8 = sc8;
    }

    @Override
    public double df(int nder, double x) {
        switch (nder) {
            //c
            case 0:
                return ((sa1*(Math.cos((L*sb1/2)-sc1)-Math.cos(sc1+sb1*x))/sb1)+(sa2*(Math.cos((L*sb2/2)-sc2)-Math.cos(sc2+sb2*x))/sb2)+(sa3*(Math.cos((L*sb3/2)-sc3)-Math.cos(sc3+sb3*x))/sb3)+(sa4*(Math.cos((L*sb4/2)-sc4)-Math.cos(sc4+sb4*x))/sb4)+(sa5*(Math.cos((L*sb5/2)-sc5)-Math.cos(sc5+sb5*x))/sb5)+(sa6*(Math.cos((L*sb6/2)-sc6)-Math.cos(sc6+sb6*x))/sb6)+(sa7*(Math.cos((L*sb7/2)-sc7)-Math.cos(sc7+sb7*x))/sb7)+(sa8*(Math.cos((L*sb8/2)-sc8)-Math.cos(sc8+sb8*x))/sb8));
            // p
            case 1:
                return (sa1*Math.sin(sb1*x+sc1)+sa2*Math.sin(sb2*x+sc2)+sa3*Math.sin(sb3*x+sc3)+sa4*Math.sin(sb4*x+sc4)+sa5*Math.sin(sb5*x+sc5)+sa6*Math.sin(sb6*x+sc6)+sa7*Math.sin(sb7*x+sc7)+sa8*Math.sin(sb8*x+sc8));
            // dp/dz
            case 2:
                return (sa1*sb1*Math.cos(sb1*x+sc1)+sa2*sb2*Math.cos(sb2*x+sc2)+sa3*sb3*Math.cos(sb3*x+sc3)+sa4*sb4*Math.cos(sb4*x+sc4)+sa5*sb5*Math.cos(sb5*x+sc5)+sa6*sb6*Math.cos(sb6*x+sc6)+sa7*sb7*Math.cos(sb7*x+sc7)+sa8*sb8*Math.cos(sb8*x+sc8));
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }



}

package etomica.meam;

import etomica.api.*;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * EAM (Embedded Atom Method) potential
 * 
 * @author Joe Kromph
 * @author Sabry Moustafa
 */
public class P2EAM extends PotentialN implements PotentialSoft {

    protected double n2 , m2 , eps , a, a2 , C , rc12, rc22;
    protected IBoundary boundary;
    protected final IVectorMutable dr, drij, drik, drjk;
    protected IVectorMutable[] gradient;
    protected IVectorMutable[] rhograd;
    protected double[] rho;

    public P2EAM(ISpace space, double n , double  m , double  eps , double  a , double  C, double rc1, double rc2) {
        super(space);
        
        n2=0.5*n;
        m2=0.5*m;
        this.eps = eps;
        this.a = a;
        a2 = a*a;
        this.C = C;
        rc12 = rc1*rc1;
        rc22 = rc2*rc2;
        dr=space.makeVector();
        drij=space.makeVector();
        drik=space.makeVector();
        drjk=space.makeVector();
        gradient=new IVectorMutable[0];
        rhograd=new IVectorMutable[0];
    }
    
    public double getRange() {
        return Math.sqrt(rc12 > rc22 ? rc12 : rc22);
    }
    public void setCutoff(double rc1, double rc2) {
        rc12 = rc1*rc1;
        rc22 = rc2*rc2;
    }
    
    public double energy(IAtomList atoms) {
        IVector pos0=atoms.getAtom(0).getPosition();
        IVector pos1=atoms.getAtom(1).getPosition();
        dr.Ev1Mv2(pos1, pos0);
        boundary.nearestImage(dr);
        double r2=dr.squared();
        double u = 0;
        if (r2<rc12) {
            u = Math.pow(a2 / r2, n2);
        }
        if (r2<rc22) {
            double rhoi = Math.pow(a2/r2, m2);
            rho[atoms.getAtom(0).getLeafIndex()] += rhoi;
            rho[atoms.getAtom(1).getLeafIndex()] += rhoi;
        }

        return u;
    }
    
    public void setBox(IBox box) {
        boundary=box.getBoundary();
        if (rho.length != box.getLeafList().getAtomCount()) {
            rho = new double[box.getLeafList().getAtomCount()];
        }
    }

    public double virial(IAtomList atoms) {
        throw new RuntimeException("implement me");
    }

    public IVector[] gradient(IAtomList atoms) {

        IVector pos0=atoms.getAtom(0).getPosition();
        IVector pos1=atoms.getAtom(1).getPosition();
        dr.Ev1Mv2(pos1, pos0);
        boundary.nearestImage(dr);
        double r2=dr.squared();
        double u = 0;
        if (r2<rc12) {
            u = Math.pow(a2 / r2, n2);
            double dvdr =  -2*eps*n2/a*Math.pow(a2/r2 , n2) ;
            gradient[1].Ea1Tv1(-dvdr, dr);
            gradient[0].Ea1Tv1(+dvdr, dr);
        }
        else {
            gradient[0].E(0);
            gradient[1].E(0);
        }
        if (r2<rc22) {
            double rhoi = Math.pow(a2/r2, m2);
            double drhodr = -2*m2*rhoi;
            rho[atoms.getAtom(0).getLeafIndex()] += rhoi;
            rho[atoms.getAtom(1).getLeafIndex()] += rhoi;
        }

        //S: Does NOT start from 0 as it will be -SUM(j>0); see below
        for(int j=1;j<atoms.getAtomCount();j++){
            gradient[j].E(0);
            rhograd[j].E(0);
            IVector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            if(rij<=rC1 && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
	            dvdr =  -eps*n/a*Math.pow(a/rij , n+1) ;
	            gradient[j].Ea1Tv1(-dvdr/rij, dr);//fji <= -fij
            }
            if(rij<=rC2){
                double drhodr;
                rhoi += Math.pow(a/rij , m);
                drhodr = -2*m2/a*Math.pow(a2/r2, m);
                rhograd[j].Ea1Tv1(-drhodr, dr);//fji <= -fij
            }          
        }
        double f=Math.sqrt(rhoi);
        for (int j=1;j<atoms.getAtomCount();j++){
            gradient[j].PEa1Tv1(-eps*C/2.0/f,rhograd[j]);//Adds the n-body to the 2-body for j
            gradient[0].ME(gradient[j]);
        }
        return gradient;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
    
}

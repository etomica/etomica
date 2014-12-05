/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.ElectronVolt;

/**
 * EFS (Extended Finnis-Sinclair) potential
 * 
 * @author Joe Kromph
 */
public class PotentialEFS extends PotentialN implements PotentialSoft{

    public static void main(String[]args){
        ISpace space=Space3D.getInstance();
        Simulation sim=new Simulation(space);
        Box box=new Box(space);
        sim.addBox(box);
        Species species=new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);
        
        box.setNMolecules(species, 2);
        
        double c,c0,c1,c2,c3,c4,A,d,B ;
        c=3.25;
        c0=ElectronVolt.UNIT.toSim(48.52796);
        c1=ElectronVolt.UNIT.toSim(-33.79621);
        c2=ElectronVolt.UNIT.toSim(5.854334);
        c3=ElectronVolt.UNIT.toSim(-0.0098221);
        c4=ElectronVolt.UNIT.toSim(0.033338);
        A=ElectronVolt.UNIT.toSim(1.885948);
        d=4.41;
        B=0;
        int qmax=1000;
        
        PotentialEFS potential=new PotentialEFS(space,c,c0,c1,c2,c3,c4,A,d,B);
        potential.setBox(box);
        
        for(int q=500; q<qmax; q++){
            box.getLeafList().getAtom(1).getPosition().setX(0, q*.01);
            double energy=potential.energy(box.getLeafList());
            System.out.println(q*.01+"   "+energy);
        } 
    }

    protected final double c, c0, c1, c2, c3, c4;
    protected final double A, d, B;
    protected IBoundary boundary; 
    protected final IVectorMutable dr;
    protected IVectorMutable[] gradient; 
    protected IVectorMutable[] rhograd;
    protected double [][] secondder;
    
    public PotentialEFS(ISpace space, double c, double c0, double c1, 
            double c2, double c3, double c4, double A, double d, double B) {
        super(space);
        
        this.c=c;
        this.c0=c0;
        this.c1=c1;
        this.c2=c2;
        this.c3=c3;
        this.c4=c4;
        this.A=A;
        this.d=d;
        this.B=B;    
        dr=space.makeVector();
        gradient=new IVectorMutable[0];
        rhograd=new IVectorMutable[0];
    }
    
    public double getRange() {
        return d;
    }
    
    public double energy(IAtomList atoms) {
        double sumV=0;
        double rhoi=0;
        IVector ipos=atoms.getAtom(0).getPosition();
        
        for(int j=1;j<atoms.getAtomCount();j++){
            IVector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            if(rij<=c && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
                sumV+=(rij-c)*(rij-c)*(c0+c1*rij+c2*rij*rij+c3*rij*rij*rij+c4*rij*rij*rij*rij);
            }
            if(rij<=d){
                double rd=rij-d;
                rhoi+=rd*rd*(1+B*B*rd*rd)*A*A;
            }
        }
        double frho=Math.sqrt(rhoi);
        double Utot=sumV-frho;
        return Utot;
        
    }
    
    public void setBox(IBox box) {
        boundary=box.getBoundary();
    }

    public double virial(IAtomList atoms) {
        double virial=0;
        gradient(atoms);
        IVector ipos=atoms.getAtom(0).getPosition();
        for(int j=1;j<atoms.getAtomCount();j++){
            IVector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            virial -= gradient[j].dot(dr);
        }
        return virial;
    }

    public IVector[] gradient(IAtomList atoms) {
        
        if(gradient.length<atoms.getAtomCount()){
            rhograd=new IVectorMutable[atoms.getAtomCount()];
            gradient=new IVectorMutable[atoms.getAtomCount()];
            
            for(int j=0;j<atoms.getAtomCount();j++){
                gradient[j]=space.makeVector();
                rhograd[j]=space.makeVector();
            }
        }
        
        gradient[0].E(0);
        IVector ipos=atoms.getAtom(0).getPosition();
        
        double rhoi=0;
        for(int j=1;j<atoms.getAtomCount();j++){
            gradient[j].E(0);
            rhograd[j].E(0);
            IVector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            if(rij<=c && atoms.getAtom(0).getLeafIndex() < atoms.getAtom(j).getLeafIndex()){
                double dvdr=(rij-c)*(2*(c0+rij*(c1+rij*(c2+rij*(c3+c4*rij))))+(rij-c)*(c1+rij*(2*c2+rij*(3*c3+4*c4*rij))));
                gradient[j].Ea1Tv1(-dvdr/rij, dr);
            }
            if(rij<=d){
                double rd=rij-d;
                rhoi+=rd*rd*(1+B*B*rd*rd)*A*A;              
                double drhodr=A*A*(2*rd+4*B*B*rd*rd*rd);
                rhograd[j].Ea1Tv1(-drhodr/rij, dr);
            }          
        }
        double f=Math.sqrt(rhoi);
        for (int j=1;j<atoms.getAtomCount();j++){
            gradient[j].PEa1Tv1(-.5*1/f,rhograd[j]);
            gradient[0].ME(gradient[j]);
        }
        return gradient;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }
    
    public double [][] secondder(IAtomList atoms){
        IVector ipos=atoms.getAtom(0).getPosition();
        
        for(int j=1;j<atoms.getAtomCount();j++){
            IVector jpos=atoms.getAtom(j).getPosition();
            dr.Ev1Mv2(ipos, jpos);
            boundary.nearestImage(dr);
            double rij=Math.sqrt(dr.squared());
            if(rij<=c){ //calculating dv^2/dxdx, dv^2/dxdx, dv^2/dxdz etc for all i or j components
                double dv2dr2=(3*rij*(2*c4*rij+c3)+c2)*(c-rij)*(c-rij)+rij*(rij*(rij*      
                        (c4*rij+c3)+c2)+c1)+2*(rij-c)*(rij*(rij*(4*c4*rij+3*c3)+2*c2)+c1)+c0;
                double dvdr=(rij-c)*(2*(c0+rij*(c1+rij*(c2+rij*(c3+c4*rij))))+(rij-c)*(c1+rij*(2*c2+rij*(3*c3+4*c4*rij))));
                for (int q=0; q<3; q++){
                    double Q= dr.getX(q);
                    double ddqqq=-(rij*rij-Q*Q)/(rij*rij*rij);
                    secondder[3*j+q][q]=ddqqq*dvdr;
                    secondder[q][q]-=secondder[3*j+q][q];
                    secondder[q][3*j+q]=secondder[3*j+q][q];
                    secondder[3*j+q][3*j+q]=-secondder[3*j+q][q];
                    for (int w=q+1; w<3; w++){
                        double W= dr.getX(w);
                        double drdqdw=-(Q*W/(rij*rij)); //we have finished filling in for components ij and ji, ii and jj.
                        double ddqqw=-(Q*W)/(rij*rij*rij);
                        secondder[3*j+q][w]=dv2dr2*drdqdw+ddqqw*dvdr;
                        secondder[3*j+w][q]=secondder[3*j+q][w];
                        secondder[w][3*j+q]=secondder[3*j+q][w];
                        secondder[q][3*j+w]=secondder[3*j+q][w];
                        secondder[w][q]-=secondder[3*j+q][w];
                        secondder[q][w]-=secondder[3*j+q][w];
                        secondder[3*j+q][3*j+w]=-secondder[3*j+q][w];
                        secondder[3*j+w][3*j+q]=-secondder[3*j+q][w];
                    }
                }
            }
        }
        return secondder;
    }
}

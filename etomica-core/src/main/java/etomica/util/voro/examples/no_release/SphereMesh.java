// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

package etomica.util.voro.examples.no_release;

import etomica.util.collections.IntArrayList;
import etomica.util.voro.*;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class SphereMesh {

    public static class WallShell implements Wall {

        public WallShell(double xc_,double yc_,double zc_,double rc,double sc) {
            this(xc_, yc_, zc_, rc, sc, -99);
        }

        public WallShell(double xc_,double yc_,double zc_,double rc,double sc,int w_id_) {
            w_id = w_id_;
            xc = xc_;
            yc = yc_;
            zc = zc_;
            lc = rc-sc;
            uc = rc+sc;
        }

        public boolean point_inside(double x,double y,double z) {
            double rsq=(x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc);
            return rsq>lc*lc&&rsq<uc*uc;
        }

        public boolean cut_cell_base(VoronoiCellBase c, double x, double y, double z) {
            double xd=x-xc,yd=y-yc,zd=z-zc,dq=xd*xd+yd*yd+zd*zd,dq2;
            if (dq>1e-5) {
                dq2=2*(Math.sqrt(dq)*lc-dq);
                dq=2*(Math.sqrt(dq)*uc-dq);
                return c.nplane(xd,yd,zd,dq,w_id)&&c.nplane(-xd,-yd,-zd,-dq2,w_id);
            }
            return true;
        }
        
        public boolean cut_cell(VoronoiCell c, double x, double y, double z) {return cut_cell_base(c,x,y,z);}
        public boolean cut_cell(VoronoiCellNeighbor c, double x, double y, double z) {return cut_cell_base(c,x,y,z);}

        private final int w_id;
		private final double xc,yc,zc,lc,uc;
    }


    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double boxl=1.2;

        // Set up the number of blocks that the container is divided into
        final int bl=14;

        // Set the number of particles that are going to be randomly introduced
        final int particles=2000;

        final int nface=11;

        int j,k,l,ll,o;
        int[] faces = new int[nface];
        double[] p = new double[3*particles];

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        Container con = new Container(-boxl,boxl,-boxl,boxl,-boxl,boxl,bl,bl,bl,false,false,false,8);

        WallShell ws = new WallShell(0,0,0,1,0.00001);
        con.add_wall(ws);

        // Randomly add particles into the container
        for (int i=0; i<particles; ) {
            double x=boxl*(2*Math.random()-1);
            double y=boxl*(2*Math.random()-1);
            double z=boxl*(2*Math.random()-1);
            double r=x*x+y*y+z*z;
            if(r>1e-5) {
                r=1/Math.sqrt(r);x*=r;y*=r;z*=r;
                con.put(i,x,y,z);
                i++;
            }
        }

        for(l=0;l<100;l++) {
            CLoopAll vl = new CLoopAll(con);
            VoronoiCell c = new VoronoiCell();
            Arrays.fill(faces, 0);
            if(vl.start()) do if(con.compute_cell(c,vl)) {
                double[] x = new double[1];
                double[] y = new double[1];
                double[] z = new double[1];
                double[] r = new double[1];
                int[] pid = new int[1];
                vl.pos(pid,x,y,z,r);
                double[] dxyz = c.centroid();
                p[3*pid[0]]=x[0]+dxyz[0];
                p[3*pid[0]+1]=y[0]+dxyz[1];
                p[3*pid[0]+2]=z[0]+dxyz[2];

                int i=c.number_of_faces()-4;
                if(i<0) i=0;if(i>=nface) i=nface-1;
                faces[i]++;
            } while (vl.inc());
            con.clear();
            double fac=0;//l<9000?0.1/sqrt(double(l)):0;
            for(int i=0;i<particles;i++) {
                con.put(i,p[3*i]+fac*(2*Math.random()-1),p[3*i+1]+fac*(2*Math.random()-1),p[3*i+2]+fac*(2*Math.random()-1));
            }
            System.out.printf("%d",l);
            for(int fp=0; fp<nface; fp++) System.out.printf(" %d",faces[fp]);
            System.out.println();
        }

        // Output the particle positions in gnuplot format
        con.draw_particles("sphere_mesh_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("sphere_mesh_v.gnu");

        // Allocate memory for neighbor relations
        int[] q=new int[particles*nface];
        int[] qn=new int[particles];
        int qp;

        // Create a table of all neighbor relations
        IntArrayList vi = new IntArrayList();
        VoronoiCellNeighbor c = new VoronoiCellNeighbor();
        CLoopAll vl = new CLoopAll(con);
        if(vl.start()) do if(con.compute_cell(c,vl)) {
            int i=vl.pid();
            qp=i*nface; // q
            c.neighbors(vi);
            if(vi.size()>nface+2) Common.voro_fatal_error("Too many faces; boost nface", Config.Voropp.INTERNAL_ERROR);

            for(l=0;l<vi.size();l++) if(vi.getInt(l)>=0) q[qp+qn[i]++]=vi.getInt(l);
        } while (vl.inc());

        // Sort the connections in anti-clockwise order
        boolean connect;
        int tote=0;
        for(l=0;l<particles;l++) {
            tote+=qn[l];
            for(int i=0;i<qn[l]-2;i++) {
                o=q[l*nface+i];
                //printf("---> %d,%d\n",i,o);
                j=i+1;
                while(j<qn[l]-1) {
                    ll=q[l*nface+j];
                    //	printf("-> %d %d\n",j,ll);
                    connect=false;
                    for(k=0;k<qn[ll];k++) {
                        //		printf("%d %d %d\n",ll,k,q[ll*nface+k]);
                        if(q[ll*nface+k]==o) {connect=true;break;}
                    }
                    if(connect) break;
                    j++;
                }

                // Swap the connected vertex into this location
                //printf("%d %d\n",i+1,j);
                o=q[l*nface+i+1];
                q[l*nface+i+1]=q[l*nface+j];
                q[l*nface+j]=o;
            }

            // Reverse all connections if the have the wrong handedness
            j=3*l;k=3*q[l*nface];o=3*q[l*nface+1];
            double x=p[j]-p[k];
            double dx=p[j]-p[o];
            double y=p[j+1]-p[k+1];
            double dy=p[j+1]-p[o+1];
            double z=p[j+2]-p[k+2];
            double dz=p[j+2]-p[o+2];
            if(p[j]*(y*dz-z*dy)+p[j+1]*(z*dx-x*dz)+p[j+2]*(x*dy-y*dx)<0) {
                for(int i=0;i<qn[l]/2;i++) {
                    o=q[l*nface+i];
                    q[l*nface+i]=q[l*nface+qn[l]-1-i];
                    q[l*nface+qn[l]-1-i]=o;
                }
            }
        }

        int[] mp = new int[particles], mpi = new int[particles];
        try {
            FileOutputStream ff = new FileOutputStream("sphere_mesh.net");
            for (int i = 1; i < particles; i++) mp[i] = -1;
            l = 1;
            o = 0;
            while (o < l) {
                int i = mpi[o];
                for (j = 0; j < qn[i]; j++) {
                    k = q[i * nface + j];
                    if (mp[k] == -1) {
                        mpi[l] = k;
                        mp[k] = l++;
                    }
                    if (mp[i] < mp[k])
                        ff.write(String.format("%g %g %g\n%g %g %g\n\n\n", p[3 * i], p[3 * i + 1], p[3 * i + 2], p[3 * k], p[3 * k + 1], p[3 * k + 2]).getBytes());
                }
                o++;
            }
            ff.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

        try {
            // Save binary representation of the mesh
            DataOutputStream fb = new DataOutputStream(new FileOutputStream("sphere_mesh.bin"));

            // Write header
            int sz = tote + particles + 2;
            fb.writeInt(1);
            fb.writeInt(sz);
            fb.writeInt(3*particles);

            // Assemble the connections and write them
            fb.writeInt(particles);
            fb.writeInt(tote);
            for (l = 0; l < particles; l++) {
                fb.writeInt(qn[mpi[l]]);
             }
            for (l = 0; l < particles; l++) {
                int i = mpi[l];
                System.out.printf("%d", l);
                for (j = 0; j < qn[i]; j++) {
                    int w = mp[q[i * nface + j]];
                    fb.writeInt(w);
                    System.out.printf(" %d", w);
                }
                System.out.println();
            }

            for (int i = 0; i < particles; i++) {
                int b = 3 * mpi[i];
                fb.writeDouble(p[b]);
                fb.writeDouble(p[b+1]);
                fb.writeDouble(p[b+2]);
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
}

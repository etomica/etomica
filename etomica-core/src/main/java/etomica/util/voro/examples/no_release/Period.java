/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples.no_release;

import etomica.util.voro.ContainerPeriodic;

public class Period {

    public static void main(String[] args) {
        // Set up constants for the container geometry
        final double bx=10;
        final double by=10;
        final double bz=10;
        final double bxy=0;
        final double bxz=5;
        final double byz=0;
        
        // Set up the number of blocks that the container is divided
        // into
        final int n_x=3,n_y=3,n_z=3;
        
        // Set the number of particles to add
        final int particles=20;

        int i;
        double x,y,z;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block.	
        ContainerPeriodic con = new ContainerPeriodic(bx,bxy,by,bxz,byz,bz,n_x,n_y,n_z,8);

        // Add particles into the container at random positions
        for(i=0;i<particles;i++) {
            x=bx*Math.random();
            y=by*Math.random();
            z=bz*Math.random();
            con.put(i,x,y,z);
        }


        // Output volume
        double vvol=con.sum_cell_volumes();
        System.out.printf("Container volume : %g\n"+
                "Voronoi volume   : %g\n",bx*by*bz,vvol);

        // Output particle positions, Voronoi cells, and the domain
        con.draw_particles("particles_periodic.gnu");
        con.draw_cells_gnuplot("cells_periodic.gnu");
        con.draw_domain_gnuplot("domain_periodic.gnu");

    }
}

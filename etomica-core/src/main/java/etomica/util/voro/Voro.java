/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

public class Voro {

    enum blocks_mode {
        none,
        length_scale,
        specified
    };

    // A maximum allowed number of regions, to prevent enormous amounts of memory
    // being allocated
    static final int max_regions=16777216;

    // This message gets displayed if the user requests the help flag
    public static void help_message() {
        System.out.println("Voro++ version 0.4.6, by Chris H. Rycroft (UC Berkeley/LBL)\n\n"+
                "Syntax: voro++ [options] <x_min> <x_max> <y_min>\n"+
                "               <y_max> <z_min> <z_max> <filename>\n\n"+
                "By default, the utility reads in the input file of particle IDs and positions,\n"+
                "computes the Voronoi cell for each, and then creates <filename.vol> with an\n"+
                "additional column containing the volume of each Voronoi cell.\n\n"+
                "Available options:\n"+
                " -c <str>   : Specify a custom output string\n"+
                " -g         : Turn on the gnuplot output to <filename.gnu>\n"+
                " -h/--help  : Print this information\n"+
                " -hc        : Print information about custom output\n"+
                " -l <len>   : Manually specify a length scale to configure the internal\n"+
                "              computational grid\n"+
                " -m <mem>   : Manually choose the memory allocation per grid block\n"+
                "              (default 8)\n"+
                " -n [3]     : Manually specify the internal grid size\n"+
                " -o         : Ensure that the output file has the same order as the input\n"+
                "              file\n"+
                " -p         : Make container periodic in all three directions\n"+
                " -px        : Make container periodic in the x direction\n"+
                " -py        : Make container periodic in the y direction\n"+
                " -pz        : Make container periodic in the z direction\n"+
                " -r         : Assume the input file has an extra coordinate for radii\n"+
                " -v         : Verbose output\n"+
                " --version  : Print version information\n"+
                " -wb [6]    : Add six plane wall objects to make rectangular box containing\n"+
                "              the space x1<x<x2, x3<y<x4, x5<z<x6\n"+
                " -wc [7]    : Add a cylinder wall object, centered on (x1,x2,x3),\n"+
                "              pointing in (x4,x5,x6), radius x7\n"+
                " -wo [7]    : Add a conical wall object, apex at (x1,x2,x3), axis\n"+
                "              along (x4,x5,x6), angle x7 in radians\n"+
                " -ws [4]    : Add a sphere wall object, centered on (x1,x2,x3),\n"+
                "              with radius x4\n"+
                " -wp [4]    : Add a plane wall object, with normal (x1,x2,x3),\n"+
                "              and displacement x4\n"+
                " -y         : Save POV-Ray particles to <filename_p.pov> and POV-Ray Voronoi\n"+
                "              cells to <filename_v.pov>\n"+
                " -yp        : Save only POV-Ray particles to <filename_p.pov>\n"+
                " -yv        : Save only POV-Ray Voronoi cells to <filename_v.pov>");
    }

    // This message gets displayed if the user requests information about doing
// custom output
    public static void custom_output_message() {
        System.out.println("The \"-c\" option allows a string to be specified that will customize the output\n"+
                "file to contain a variety of statistics about each computed Voronoi cell. The\n"+
                "string is similar to the standard C printf() function, made up of text with\n"+
                "additional control sequences that begin with percentage signs that are expanded\n"+
                "to different statistics. See http://math.lbl.gov/voro++/doc/custom.html for more\n"+
                "information.\n"+
                "\nParticle-related:\n"+
                "  %i The particle ID number\n"+
                "  %x The x coordinate of the particle\n"+
                "  %y The y coordinate of the particle\n"+
                "  %z The z coordinate of the particle\n"+
                "  %q The position vector of the particle, short for \"%x %y %z\"\n"+
                "  %r The radius of the particle (only printed if -p enabled)\n"+
                "\nVertex-related:\n"+
                "  %w The number of vertices in the Voronoi cell\n"+
                "  %p A list of the vertices of the Voronoi cell in the format (x,y,z),\n"+
                "     relative to the particle center\n"+
                "  %P A list of the vertices of the Voronoi cell in the format (x,y,z),\n"+
                "     relative to the global coordinate system\n"+
                "  %o A list of the orders of each vertex\n"+
                "  %m The maximum radius squared of a vertex position, relative to the\n"+
                "     particle center\n"+
                "\nEdge-related:\n"+
                "  %g The number of edges of the Voronoi cell\n"+
                "  %E The total edge distance\n"+
                "  %e A list of perimeters of each face\n"+
                "\nFace-related:\n"+
                "  %s The number of faces of the Voronoi cell\n"+
                "  %F The total surface area of the Voronoi cell\n"+
                "  %A A frequency table of the number of edges for each face\n"+
                "  %a A list of the number of edges for each face\n"+
                "  %f A list of areas of each face\n"+
                "  %t A list of bracketed sequences of vertices that make up each face\n"+
                "  %l A list of normal vectors for each face\n"+
                "  %n A list of neighboring particle or wall IDs corresponding to each face\n"+
                "\nVolume-related:\n"+
                "  %v The volume of the Voronoi cell\n"+
                "  %c The centroid of the Voronoi cell, relative to the particle center\n"+
                "  %C The centroid of the Voronoi cell, in the global coordinate system");
    }


    // Ths message is displayed if the user requests version information
    public static void version_message() {
        System.out.println("Voro++ version 0.4.6 (October 17th 2013)");
    }

    // Prints an error message. This is called when the program is unable to make
    // sense of the command-line options.
    public static void error_message() {
        System.err.println("voro++: Unrecognized command-line options; type \"voro++ -h\" for more\ninformation.");
    }

    // Carries out the Voronoi computation and outputs the results to the requested
// files
    public static void cmd_line_output(CLoopBase vl, ContainerBase con, String format, OutputStream outfile, OutputStream gnu_file, OutputStream povp_file, OutputStream povv_file, boolean verbose, double[] vol, int[] vcc, int[] tp) {
        int[] pid = new int[1];
        int ps=con.ps;
        double[] x = new double[1];
        double[] y = new double[1];
        double[] z = new double[1];
        double[] r = new double[1];
        if(con.contains_neighbor(format)) {
            if(vl.start()) {
                VoronoiCellNeighbor c = new VoronoiCellNeighbor(con);
                VoronoiCellNeighbor[] cout = new VoronoiCellNeighbor[]{c};
                do {
                    if(con.compute_cell(cout,vl)) {
                        vl.pos(pid,x,y,z,r);
                        if(outfile!=null) cout[0].output_custom(format,pid[0],x[0],y[0],z[0],r[0],outfile);
                        if(gnu_file!=null) cout[0].draw_gnuplot(x[0],y[0],z[0],gnu_file);
                        if(povp_file!=null) {
                            try {
                                povp_file.write(String.format("// id %d\n", pid[0]).getBytes());
                                if (ps == 4)
                                    povp_file.write(String.format("sphere{<%g,%g,%g>,%g}\n", x[0], y[0], z[0], r[0]).getBytes());
                                else
                                    povp_file.write(String.format("sphere{<%g,%g,%g>,s}\n", x[0], y[0], z[0]).getBytes());
                            }
                            catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                        }
                        if(povv_file!=null) {
                            try {
                                povv_file.write(String.format("// cell %d\n", pid[0]).getBytes());
                            }
                            catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                            cout[0].draw_pov(x[0],y[0],z[0],povv_file);
                        }
                        if(verbose) {
                            vol[0]+=c.volume();
                            vcc[0]++;
                        }
                    }
                } while(vl.inc());
            }
        } else {
            if(vl.start()) {
                VoronoiCell c = new VoronoiCell(con);
                VoronoiCell[] cout = new VoronoiCell[]{c};
                do {
                    if(con.compute_cell(cout,vl)) {
                        vl.pos(pid,x,y,z,r);
                        if(outfile!=null) cout[0].output_custom(format,pid[0],x[0],y[0],z[0],r[0],outfile);
                        if(gnu_file!=null) cout[0].draw_gnuplot(x[0],y[0],z[0],gnu_file);
                        if(povp_file!=null) {
                            try {
                                povp_file.write(String.format("// id %d\n", pid[0]).getBytes());
                                if (ps == 4)
                                    povp_file.write(String.format("sphere{<%g,%g,%g>,%g}\n", x[0], y[0], z[0], r[0]).getBytes());
                                else
                                    povp_file.write(String.format("sphere{<%g,%g,%g>,s}\n", x[0], y[0], z[0]).getBytes());
                            }
                            catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                        }
                        if(povv_file!=null) {
                            try {
                                povv_file.write(String.format("// cell %d\n", pid[0]).getBytes());
                                cout[0].draw_pov(x[0], y[0], z[0], povv_file);
                            }
                            catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                        }
                        if(verbose) {
                            vol[0]+=cout[0].volume();
                            vcc[0]++;
                        }
                    }
                } while(vl.inc());
            }
        }
        if(verbose) tp[0]=con.total_particles();
    }

    public static void main(String[] args) {
        int j=-7,custom_output=0;
        int nx=0,ny=0,nz=0,init_mem=8;
        double ls=0;
        blocks_mode bm=blocks_mode.none;
        boolean gnuplot_output=false,povp_output=false,povv_output=false,polydisperse=false;
        boolean xperiodic=false,yperiodic=false,zperiodic=false,ordered=false,verbose=false;
        PreContainer pcon = null;
        PreContainerPoly pconp = null;
        ArrayList<Wall> wl = new ArrayList<>();

        // If there's one argument, check to see if it's requesting help.
        // Otherwise, bail out with an error.
        if(args.length==1) {
            if(args[0].equals("-h") || args[0].equals("--help")) {
                help_message();
                return;
            } else if(args[0].equals("-hc")) {
                custom_output_message();
                return;
            } else if(args[0].equals("--version")) {
                version_message();
                return;
            } else {
                error_message();
                throw new RuntimeException();
            }
        }

        // If there aren't enough command-line arguments, then bail out
        // with an error.
        if(args.length<6) {
            error_message();
            throw new RuntimeException();
        }

        // We have enough arguments. Now start searching for command-line
        // options.
        int i;
        for (i=0; i<args.length-6; i++) {
            if(args[i].equals("-c")) {
                if(i>=args.length-7) {
                    error_message();
                    throw new RuntimeException();
                }
                if(custom_output==0) {
                    custom_output=++i;
                } else {
                    throw new RuntimeException("voro++: multiple custom output strings detected\n");
                }
            } else if(args[i].equals("-g")) {
                gnuplot_output=true;
            } else if(args[i].equals("-h") || args[i].equals("--help")) {
                help_message();
                return;
            } else if(args[i].equals("-hc")) {
                custom_output_message();
                return;
            } else if(args[i].equals("-l")) {
                if(i>=args.length-7) {
                    error_message();
                    throw new RuntimeException();
                }
                if(bm!=blocks_mode.none) {
                    throw new RuntimeException("voro++: Conflicting options about grid setup (-l/-n)");
                }
                bm=blocks_mode.length_scale;
                i++;
                ls=Double.parseDouble(args[i]);
            } else if(args[i].equals("-m")) {
                i++;
                init_mem=Integer.parseInt(args[i]);
            } else if(args[i].equals("-n")) {
                if(i>=args.length-9) {
                    error_message();
                    throw new RuntimeException();
                }
                if(bm!=blocks_mode.none) {
                    throw new RuntimeException("voro++: Conflicting options about grid setup (-l/-n)");
                }
                bm=blocks_mode.specified;
                i++;
                nx=Integer.parseInt(args[i++]);
                ny=Integer.parseInt(args[i++]);
                nz=Integer.parseInt(args[i]);
                if(nx<=0||ny<=0||nz<=0) {
                    throw new RuntimeException("voro++: Computational grid specified with -n must be greater than one\n"+
                            "in each direction");
                }
            } else if(args[i].equals("-o")) {
                ordered=true;
            } else if(args[i].equals("-p")) {
                xperiodic=yperiodic=zperiodic=true;
            } else if(args[i].equals("-px")) {
                xperiodic=true;
            } else if(args[i].equals("-py")) {
                yperiodic=true;
            } else if(args[i].equals("-pz")) {
                zperiodic=true;
            } else if(args[i].equals("-r")) {
                polydisperse=true;
            } else if(args[i].equals("-v")) {
                verbose=true;
            } else if(args[i].equals("--version")) {
                version_message();
                throw new RuntimeException();
            } else if(args[i].equals("-wb")) {
                if(i>=args.length-12) {
                    error_message();
                    throw new RuntimeException();
                }
                i++;
                double w0=Double.parseDouble(args[i++]);
                double w1=Double.parseDouble(args[i++]);
                double w2=Double.parseDouble(args[i++]);
                double w3=Double.parseDouble(args[i++]);
                double w4=Double.parseDouble(args[i++]);
                double w5=Double.parseDouble(args[i]);
                wl.add(new WallPlane(-1,0,0,-w0,j));j--;
                wl.add(new WallPlane(1,0,0,w1,j));j--;
                wl.add(new WallPlane(0,-1,0,-w2,j));j--;
                wl.add(new WallPlane(0,1,0,w3,j));j--;
                wl.add(new WallPlane(0,0,-1,-w4,j));j--;
                wl.add(new WallPlane(0,0,1,w5,j));j--;
            } else if(args[i].equals("-ws")) {
                if(i>=args.length-10) {
                    error_message();
                    throw new RuntimeException();
                }
                i++;
                double w0=Double.parseDouble(args[i++]);
                double w1=Double.parseDouble(args[i++]);
                double w2=Double.parseDouble(args[i++]);
                double w3=Double.parseDouble(args[i]);
                wl.add(new WallSphere(w0,w1,w2,w3,j));
                j--;
            } else if(args[i].equals("-wp")) {
                if(i>=args.length-10) {
                    error_message();
                    throw new RuntimeException();
                }
                i++;
                double w0 = Double.parseDouble(args[i++]);
                double w1 = Double.parseDouble(args[i++]);
                double w2 = Double.parseDouble(args[i++]);
                double w3 = Double.parseDouble(args[i]);
                wl.add(new WallPlane(w0,w1,w2,w3,j));
                j--;
            } else if(args[i].equals("-wc")) {
                if(i>=args.length-13) {
                    error_message();
                    throw new RuntimeException();
                }
                i++;
                double w0 = Double.parseDouble(args[i++]);
                double w1 = Double.parseDouble(args[i++]);
                double w2 = Double.parseDouble(args[i++]);
                double w3 = Double.parseDouble(args[i++]);
                double w4 = Double.parseDouble(args[i++]);
                double w5 = Double.parseDouble(args[i++]);
                double w6 = Double.parseDouble(args[i]);
                wl.add(new WallCylinder(w0,w1,w2,w3,w4,w5,w6,j));
                j--;
            } else if(args[i].equals("-wo")) {
                if(i>=args.length-13) {
                    error_message();
                    throw new RuntimeException();
                }
                i++;
                double w0 = Double.parseDouble(args[i++]);
                double w1 = Double.parseDouble(args[i++]);
                double w2 = Double.parseDouble(args[i++]);
                double w3 = Double.parseDouble(args[i++]);
                double w4 = Double.parseDouble(args[i++]);
                double w5 = Double.parseDouble(args[i++]);
                double w6 = Double.parseDouble(args[i]);
                wl.add(new WallCone(w0,w1,w2,w3,w4,w5,w6,j));
                j--;
            } else if(args[i].equals("-y")) {
                povp_output=povv_output=true;
            } else if(args[i].equals("-yp")) {
                povp_output=true;
            } else if(args[i].equals("-yv")) {
                povv_output=true;
            } else {
                error_message();
                throw new RuntimeException();
            }
            i++;
        }

        // Check the memory guess is positive
        if(init_mem<=0) {
            throw new RuntimeException("voro++: The memory allocation must be positive");
        }

        // Read in the dimensions of the test box, and estimate the number of
        // boxes to divide the region up into
        double ax=Double.parseDouble(args[i]);
        double bx=Double.parseDouble(args[i+1]);
        double ay=Double.parseDouble(args[i+2]);
        double by=Double.parseDouble(args[i+3]);
        double az=Double.parseDouble(args[i+4]);
        double bz=Double.parseDouble(args[i+5]);

        // Check that for each coordinate, the minimum value is smaller
        // than the maximum value
        if(bx<ax) {
            throw new RuntimeException("voro++: Minimum x coordinate exceeds maximum x coordinate");
        }
        if(by<ay) {
            throw new RuntimeException("voro++: Minimum y coordinate exceeds maximum y coordinate");
        }
        if(bz<az) {
            throw new RuntimeException("voro++: Minimum z coordinate exceeds maximum z coordinate");
        }

        if(bm==blocks_mode.none) {
            if(polydisperse) {
                pconp=new PreContainerPoly(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
                pconp.import_(args[i+6]);
                int[] nxout = new int[]{nx};
                int[] nyout = new int[]{ny};
                int[] nzout = new int[]{nz};
                pconp.guess_optimal(nxout,nyout,nzout);
                nx = nxout[0];
                ny = nyout[0];
                nz = nzout[0];
            } else {
                pcon=new PreContainer(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
                pcon.import_(args[i+6]);
                int[] nxout = new int[]{nx};
                int[] nyout = new int[]{ny};
                int[] nzout = new int[]{nz};
                pcon.guess_optimal(nxout,nyout,nzout);
                nx = nxout[0];
                ny = nyout[0];
                nz = nzout[0];
            }
        } else {
            double nxf,nyf,nzf;
            if(bm==blocks_mode.length_scale) {

                // Check that the length scale is positive and
                // reasonably large
                if(ls<Config.tolerance) {
                    String msg = "voro++: ";
                    if(ls<0) {
                        msg += "The length scale must be positive";
                    } else {
                        msg += String.format("The length scale is smaller than the safe limit of %g. Either\nincrease the particle length scale, or recompile with a different limit.", Config.tolerance);
                    }
                    throw new RuntimeException(msg);
                }
                ls=0.6/ls;
                nxf=(bx-ax)*ls+1;
                nyf=(by-ay)*ls+1;
                nzf=(bz-az)*ls+1;

                nx=(int)nxf;ny=(int)nyf;nz=(int)nzf;
            } else {
                nxf=nx;nyf=ny;nzf=nz;
            }

            // Compute the number regions based on the length scale
            // provided. If the total number exceeds a cutoff then bail
            // out, to prevent making a massive memory allocation. Do this
            // test using floating point numbers, since huge integers could
            // potentially wrap around to negative values.
            if(nxf*nyf*nzf>Voro.max_regions) {
                String msg = String.format("voro++: Number of computational blocks exceeds the maximum allowed of %d.\n"+
                        "Either increase the particle length scale, or recompile with an increased\nmaximum.",max_regions);
                throw new RuntimeException(msg);
            }
        }

        // Check that the output filename is a sensible length
        int flen=args[i+6].length();
        if(flen>4096) {
            throw new RuntimeException("voro++: Filename too long");
        }

        // Open files for output
        String buffer = String.format("%s.vol", args[i+6]);
        FileOutputStream outfile = null;
        FileOutputStream gnu_file = null, povp_file = null, povv_file = null;
        try {
            outfile = new FileOutputStream(buffer);
            if (gnuplot_output) {
                buffer = String.format("%s.gnu", args[i + 6]);
                gnu_file = new FileOutputStream(buffer);
            }
            if (povp_output) {
                buffer = String.format("%s_p.pov", args[i + 6]);
                povp_file = new FileOutputStream(buffer);
            }
            if (povv_output) {
                buffer = String.format("%s_v.pov", args[i+6]);
                povv_file = new FileOutputStream(buffer);
            }
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    	String c_str=(custom_output==0?(polydisperse?"%i %q %v %r":"%i %q %v"):args[custom_output]);

        // Now switch depending on whether polydispersity was enabled, and
        // whether output ordering is requested
        double[] vol=new double[1];
        int[] tp=new int[1],vcc=new int[1];
        if(polydisperse) {
            if(ordered) {
                ParticleOrder vo = new ParticleOrder();
                ContainerPoly con = new ContainerPoly(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);
                con.add_wall(wl);
                if(bm==blocks_mode.none) {
                    pconp.setup(vo,con);
                } else con.import_(vo,args[i+6]);

                CLoopOrder vlo = new CLoopOrder(con,vo);
                cmd_line_output(vlo,con,c_str,outfile,gnu_file,povp_file,povv_file,verbose,vol,vcc,tp);
            } else {
                ContainerPoly con = new ContainerPoly(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);
                con.add_wall(wl);

                if(bm==blocks_mode.none) {
                    pconp.setup(con);
                } else con.import_(args[i+6]);

                CLoopAll vla = new CLoopAll(con);
                cmd_line_output(vla,con,c_str,outfile,gnu_file,povp_file,povv_file,verbose,vol,vcc,tp);
            }
        } else {
            if(ordered) {
                ParticleOrder vo = new ParticleOrder();
                Container con = new Container(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);
                con.add_wall(wl);
                if(bm==blocks_mode.none) {
                    pcon.setup(vo,con);
                } else con.import_(vo,args[i+6]);

                CLoopOrder vlo = new CLoopOrder(con,vo);
                cmd_line_output(vlo,con,c_str,outfile,gnu_file,povp_file,povv_file,verbose,vol,vcc,tp);
            } else {
                Container con = new Container(ax,bx,ay,by,az,bz,nx,ny,nz,xperiodic,yperiodic,zperiodic,init_mem);
                con.add_wall(wl);
                if(bm==blocks_mode.none) {
                    pcon.setup(con);
                } else con.import_(args[i+6]);
                CLoopAll vla = new CLoopAll(con);
                cmd_line_output(vla,con,c_str,outfile,gnu_file,povp_file,povv_file,verbose,vol,vcc,tp);
            }
        }

        // Print information if verbose output requested
        if(verbose) {
            System.out.printf("Container geometry        : [%g:%g] [%g:%g] [%g:%g]\n"+
                    "Computational grid size   : %d by %d by %d (%s)\n"+
                    "Filename                  : %s\n"+
                    "Output string             : %s%s\n",ax,bx,ay,by,az,bz,nx,ny,nz,
                    bm==blocks_mode.none?"estimated from file":(bm==blocks_mode.length_scale?
                            "estimated using length scale":"directly specified"),
                    args[i+6],c_str,custom_output==0?" (default)":"");
            System.out.printf("Total imported particles  : %d (%.2g per grid block)\n"+
                    "Total V. cells computed   : %d\n"+
                    "Total container volume    : %g\n"+
                    "Total V. cell volume      : %g\n",tp[0],((double) tp[0])/(nx*ny*nz),
                    vcc[0],(bx-ax)*(by-ay)*(bz-az),vol[0]);
        }

        // Close output files
        try {
            outfile.close();
            if (gnu_file != null) gnu_file.close();
            if (povp_file != null) povp_file.close();
            if (povv_file != null) povv_file.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }
}

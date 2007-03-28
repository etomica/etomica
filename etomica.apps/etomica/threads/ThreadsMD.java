package etomica.threads;

/**
 * 
 * @author msellers
 *
 */
		public class ThreadsMD extends Thread {
			
			boolean run = true;
			int stepcounter1 = 0;
			int stepcounter2 = 0;
			final double deltat;
			final int box;
			final int threadID;
			final int Tatoms;
			final int natoms;
			final int startAtom; 
			final int stopAtom;
			final AtomDecompositionMD AD;
			
			final double[] r = new double[3];
			final double[][] Tcoord;
			final double[][] Tvelocity;
			final double[][] Tforce;
			
			final double[][] coord;
			final double[][] velocity;
			final double[][] force; 
			
			public ThreadsMD(int t, int n, int z, int b, double dt, double[][] c, double[][] v, double[][] f, AtomDecompositionMD AD){
				
				this.AD = AD;
				this.coord = c;
				this.velocity = v;
				this.force = f;
				this.threadID = t;
				this.Tatoms = n;
				this.natoms = z;
				this.box = b;
				this.deltat = dt;
				
				startAtom = threadID * Tatoms;
				stopAtom = (threadID + 1) * Tatoms;
				
				Tcoord = new double[natoms][3];
				Tvelocity = new double[natoms][3];
				Tforce = new double[natoms][3];
			}
			
			public void run() {

				while (run){
					
					// FIRST Wait on MAIN thread
					synchronized(AD){
						try{
							if 	(!run){ break;}
							
							if (stepcounter1==AD.ready1){
								AD.wait();	
							}
						}
						catch(InterruptedException e){
							if(!run){ break;}
						}
					
					}
					
					// Run Verlet integration scheme
					Verlet();
					
					synchronized(this){
						stepcounter1++;
						notifyAll();		
					}
					
					// SECOND Wait on MAIN thread
					synchronized(AD){
						try{
							if 	(!run){ break;}
							
							if (stepcounter2==AD.ready2){
								AD.wait();	
							}
						}
						catch(InterruptedException e){
							if(!run){ break;}
						}
					
					}
					
					// Write data to coord,velocity,force
					writeData();
					
					synchronized(this){
						stepcounter2++;
						notifyAll();		
					}
				
					
					
					
					
					
					
				}
				
			}
								
			public void Verlet(){
				
				//	Compute positions and half update velocities
				for (int i=startAtom; i<stopAtom; i++){
					for (int l=0; l<3; l++){
						Tcoord[i][l] = deltat*velocity[i][l] + deltat*deltat*force[i][l] * 0.5;
						Tvelocity[i][l] = velocity[i][l] + deltat*force[i][l] * 0.5;
					}
				}
				
				//  Find current r, and calculate force
				for (int i=startAtom; i<stopAtom; i++){
					for (int j=0; j<natoms; j++){
						
						if (j-i==0){continue;}
						
						double r2 = 0.0;
						
						for(int l=0; l<3; l++){
							Tforce[i][l] = 0.0;
							
							r[l] = coord[i][l] - coord[j][l];
							r[l] = r[l] - box*Math.round(r[l]/box);
							r2 += r[l]*r[l];
						}
						
						if (r2 < 6.25){
							calcForce(i,r2);
						}
						
					}
				}
				
				for (int i=startAtom; i<stopAtom; i++){
					for(int l=0; l<3; l++){
						Tvelocity[i][l] += deltat*Tforce[i][l]*0.5;
					}
				}
				
			}
			
			public void calcForce(int i, double rsq){
			
				double r2i = 1.0/rsq;
				double r6i = r2i*r2i*r2i;
				double ff = 48*r2i*r6i*(r6i-0.5);
				for(int l=0; l<3; l++) {
					Tforce[i][l] += ff * r[l];
				}
			}

			public void writeData(){
				//  Write to global coord, velocity, and force arrays
				for (int i=startAtom; i<stopAtom; i++){
					for (int l=0; l<3; l++){
						coord[i][l] += Tcoord[i][l];
						velocity[i][l] = Tvelocity[i][l];
						force[i][l] = Tforce[i][l];
					}
				}
			}
		}



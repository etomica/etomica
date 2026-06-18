/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

import java.util.Arrays;

/** \brief A class for storing ordering information when particles are added to
 * a container.
 *
 * When particles are added to a container class, they are sorted into an
 * internal computational grid of blocks. The particle_order class provides a
 * mechanism for remembering which block particles were sorted into. The import
 * and put routines in the container class have variants that also take a
 * particle_order class. Each time they are called, they will store the block
 * that the particle was sorted into, plus the position of the particle within
 * the block. The particle_order class can used by the c_loop_order class to
 * specifically loop over the particles that have their information stored
 * within it. */
public class ParticleOrder {
    /** A pointer to the array holding the ordering. */
    int[] o;
    /** A pointer to the next position in the ordering array in
     * which to store an entry. */
    int op;
    /** The particle_order constructor allocates memory to store the
     * ordering information.
     * \param[in] init_size the initial amount of memory to
     *                      allocate. */
    public ParticleOrder() {
        this(Config.init_ordering_size);
    }
    public ParticleOrder(int init_size) {
        o = new int[2*init_size];
        op = 0;
    }
    /** Adds a record to the order, corresponding to the memory
     * address of where a particle was placed into the container.
     * \param[in] ijk the block into which the particle was placed.
     * \param[in] q the position within the block where the
     * 		particle was placed. */
    public void add(int ijk,int q) {
        if(op==o.length) o = Arrays.copyOf(o, o.length*2);
        o[op] = ijk;
        o[op+1] = q;
        op+=2;
    }
}

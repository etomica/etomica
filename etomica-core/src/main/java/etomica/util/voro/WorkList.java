/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro;

public class WorkList {

    /** Each region is divided into a grid of subregions, and a worklist is
     # constructed for each. This parameter sets is set to half the number of
     # subregions that the block is divided into. */
    public static final int wl_hgrid=4;
    /** The number of subregions that a block is subdivided into, which is twice
     the value of hgrid. */
    public static final int wl_fgrid=8;
    /** The total number of worklists, set to the cube of hgrid. */
    public static final int wl_hgridcu=64;
    /** The number of elements in each worklist. */
    public static final int wl_seq_length=64;

}

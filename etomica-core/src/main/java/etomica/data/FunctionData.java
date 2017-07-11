/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

/**
 * Interface for the basic features of a general function that maps an object
 * into Data.
 */
public interface FunctionData<T> {

    /**
     * Method that performs the mapping defined by the function.
     * 
     * @param x
     *            the input object, typically a Double, Vector, Tensor, Atom,
     *            etc.
     * @return the output instance of Data defined by this function
     */
    public IData f(T x);

    /**
     * Returns a DataInfo that describes the Data object that is returned by this function's f method.
     */
    public IDataInfo getDataInfo();

}

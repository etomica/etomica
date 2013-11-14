package g3dsys.control;

public interface IndexIterator {
  
  /**
   * Get the number of dimensions
   * @return returns the number of dimensions
   */
  public int getD();
  /**
   * Resets iteration; should be called before use
   */
  public void reset();
  /**
   * See if there are more indices
   * @return returns whether there are more indices
   */
  public boolean hasNext();
  /**
   * Get the next set of indices
   * @return returns the next indices
   */
  public int[] next();
  
  /**
   * Sets the size of each dimension; in general they should all be the
   * same size.
   * @param size an array of sizes to set
   */
  public void setSize(int[] size);
  
  /**
   * Indicates whether it is appropriate to simply scale the original
   * boundary lines when in 'large box' mode.
   * @return returns whether the boundary is lazy-safe
   */
  public boolean isLazySafe();

}

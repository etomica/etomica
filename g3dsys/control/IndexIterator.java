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
  
  public void setSize(int[] size);

}

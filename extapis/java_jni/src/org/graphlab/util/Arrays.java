package org.graphlab.util;

import java.lang.reflect.Array;

/**
 * Compatibility fix between JDK 5 and JDK 6.
 * @author Jiunn Haur Lim
 * @see <a href="http://hg.openjdk.java.net/jdk6/jdk6/jdk/raw-file/ffa98eed5766/src/share/classes/java/util/Arrays.java">Open JDK</a>
 */
public class Arrays {

  // Cloning
  /**
   * Copies the specified array, truncating or padding with nulls (if necessary)
   * so the copy has the specified length.  For all indices that are
   * valid in both the original array and the copy, the two arrays will
   * contain identical values.  For any indices that are valid in the
   * copy but not the original, the copy will contain <tt>null</tt>.
   * Such indices will exist if and only if the specified length
   * is greater than that of the original array.
   * The resulting array is of exactly the same class as the original array.
   *
   * @param original the array to be copied
   * @param newLength the length of the copy to be returned
   * @return a copy of the original array, truncated or padded with nulls
   *     to obtain the specified length
   * @throws NegativeArraySizeException if <tt>newLength</tt> is negative
   * @throws NullPointerException if <tt>original</tt> is null
   * @since 1.6
   */
  @SuppressWarnings("unchecked")
  public static <T> T[] copyOf(T[] original, int newLength) {
      return (T[]) copyOf(original, newLength, original.getClass());
  }

  /**
   * Copies the specified array, truncating or padding with nulls (if necessary)
   * so the copy has the specified length.  For all indices that are
   * valid in both the original array and the copy, the two arrays will
   * contain identical values.  For any indices that are valid in the
   * copy but not the original, the copy will contain <tt>null</tt>.
   * Such indices will exist if and only if the specified length
   * is greater than that of the original array.
   * The resulting array is of the class <tt>newType</tt>.
   *
   * @param original the array to be copied
   * @param newLength the length of the copy to be returned
   * @param newType the class of the copy to be returned
   * @return a copy of the original array, truncated or padded with nulls
   *     to obtain the specified length
   * @throws NegativeArraySizeException if <tt>newLength</tt> is negative
   * @throws NullPointerException if <tt>original</tt> is null
   * @throws ArrayStoreException if an element copied from
   *     <tt>original</tt> is not of a runtime type that can be stored in
   *     an array of class <tt>newType</tt>
   * @since 1.6
   */
  @SuppressWarnings("unchecked")
  public static <T,U> T[] copyOf(U[] original, int newLength, Class<? extends T[]> newType) {
      T[] copy = ((Object)newType == (Object)Object[].class)
          ? (T[]) new Object[newLength]
          : (T[]) Array.newInstance(newType.getComponentType(), newLength);
      System.arraycopy(original, 0, copy, 0,
                       Math.min(original.length, newLength));
      return copy;
  }

}

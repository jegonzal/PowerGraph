package org.graphlab.util;

import java.util.Collection;
import java.util.Iterator;

public class StringUtils {
  
  public static String join(Collection<String> s, String delimiter) {
    if (s == null || s.isEmpty())
      return "";
    Iterator<String> iter = s.iterator();
    StringBuilder builder = new StringBuilder(iter.next());
    while (iter.hasNext()) {
      builder.append(delimiter).append(iter.next());
    }
    return builder.toString();
  }

}

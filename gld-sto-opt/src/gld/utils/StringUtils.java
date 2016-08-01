package gld.utils;

import java.util.NoSuchElementException;

/** Auxiliary class for processing strings. Used mainly by the XML parser
 */

public class StringUtils
{
  /** Replace instances of strings in the from-array in the first string
    * with their counterpart in the to-array.
    * @param string The input string
    * @param from The first column of the translation table (The strings to
    *        replace).
    * @param to The second column of the translation table (The strings to
             replace the elements of the from array with).
    * @returns The resulting string with all stuff replaced
   */
  public static String replaceList (String string,String[] from,String[] to)
  { String result=new String(string);
    for (int t=0;t<from.length;t++)
        result=StringUtils.replace(result,from[t],to[t]);
    return result;
  }

  /** Replace instances of chars in the from-array in the first string
    * with their counterpart in the to-array. This method is way more efficient
    * than the replaceList(String,String[],String[]) method.
    * @param string The input string
    * @param from The first column of the translation table (The chars to
    *        replace).
    * @param to The second column of the translation table (The strings to
             replace the elements of the from array with).
    * @returns The resulting string with all stuff replaced
   */
  public static String replaceList (String string,char[] from,String[] to )
  {	String result="";
  	char c;
	int lastCopy=0;
	for (int t=0;t<string.length();t++)
	{	c=string.charAt(t);
		for (int u=0;u<from.length;u++)
		{	if (c==from[u])
			{	result+=string.substring(lastCopy,t)+to[u];
				lastCopy=t+1;
				break;
			}
		}
	}
	if (lastCopy<string.length())
       		result+=string.substring(lastCopy,string.length());
    	return result;
  }

  /** Replace all occurences of a string in another String
    * @param string The string to search in
    * @param o the old substring
    * @param n the string with which o has to be replaced
    * @return The result
   */
  public static String replace (String string,String o,String n)
  { String result="";
    int newIndex=0,oldIndex=0;
    while ( (newIndex=string.indexOf(o,oldIndex))!=-1)
    { result+=string.substring(oldIndex,newIndex)+n;
      oldIndex=newIndex+o.length();
    }
    result+=string.substring(oldIndex,string.length());
    return result;
  }

  /** Replace all occurences of a char in a String with a string
    * @param string The string to search in
    * @param o the old char
    * @param n the new string to replace the char with
    * @return The result
   */
  public static String replace (String string,char o,String n)
  { return StringUtils.replace(string,o+"",n);
  }

  /** Replace a char in a string with another char
    * @param string The string to search in
    * @param o the old char
    * @param n the new char with which the old char is to be replaced
    * @return The result
   */
  public static String replace (String string,char o,char n)
  { return string.replace(o,n);
  }

  /** Remove all occurences of a char from a string
    * @param string The string to search in
    * @param o the char to remove
    * @return The result
   */
  public static String remove (String string,char o)
  { return StringUtils.replace (string,o+"","");
  }

  /** Remove all occurences of a string from another string
    * @param string The string to search in
    * @param o the string to remove
    * @return The result
   */
  public static String remove (String string,String o)
  { return StringUtils.replace (string,o,"");
  }

  /** Fill a string with several occurences of a char
    * @param c The char to fill the string with
    * @param times The length of the string
    * @return A string which consists of "times" times char c
   */
  public static String repeat (char c,int times)
  { String result="";
    while (times > 0 )
    { if ( (times % 2)==0 && result.length()>0)
      { result+=result;
        times/=2;
      }
      else
      { result+=c;
        times--;
      }
    }
    return result;
  }

  /** Convert the first letter of a string to upper case and return the result
    * @param string The input string
    * @return The result
   */
  public static String firstLetterToUpperCase (String string)
  { if (string.length() < 2)
    { return string.toUpperCase();
    }
    else
    { return Character.toUpperCase(string.charAt(0))+
             string.substring(1,string.length());
    }
  }

  	/** Lookup an element in an array
	  * @param data The array to search in
	  * @param index The index of the element you want
	  * @returns The element at that position
	  * @throws NoSuchElementException If the array has no element at
	  *	    the specified index.
	 */
	public static Object lookUpNumber (Object[] data,int index)
  	{ 	if (index<0 || index >= data.length )
    			throw new NoSuchElementException
			("Stringutils : number lookup in array failed for "+
			  "number "+index);
		else
			return data[index];
	}

	/** Find the index of an element in an array. Elements are compared
	  * using the object.equals(otherObject) method.
	  * @param data The array to search in
	  * @param object The element to search for
	  * @returns The index of that element
	  * @throws NullPointerException If the object to search for is null
	  * @throws NoSuchElementException If the specified object cannot be
	  *         found
	 */
	public static int getIndexObject (Object[] data,Object object)
	{	for (int x=0;x<data.length;x++)
		{	if (object.equals(data[x]))
				return x;
		}
		throw new NoSuchElementException
		("StringUtils : index lookup in array failed for object "+
		  object+". Data length is "+data.length);
	}


}


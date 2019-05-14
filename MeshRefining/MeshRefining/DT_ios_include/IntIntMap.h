/*
* File: IntIntMap.h
* -----------
* This interface exports a slightly simplified version of the Map
* class that can replace the function of std;;map<int,int> in a very efficent manner.
*/
#ifndef _int_int_map_h
#define _int_int_map_h

/*
* Class: IntIntMap
* ----------
* This interface defines a class  that stores a collection
* of key-value pairs. Both the keys and the values are always integars. 
*/

class IntIntMap 
{
public:
	/*
	* Constructor
	* The constructor initializes a new empty map.
	*/
	IntIntMap();
	IntIntMap(int initSize);

	/*
	* Destructor: ~IntIntMap
	* The destructor frees any heap storage associated with this map.
	*/
	~IntIntMap();

	/*
	* Method: size
	* Usage: nEntries = map.size();
	* -----------------------------
	* This method returns the number of entries in this map.
	*/
	int size();

	/*
	* Method: isEmpty
	* Usage: if (map.isEmpty())...
	* ----------------------------
	* This method returns true if this map contains no entries,
	* false otherwise.
	*/
	bool isEmpty();

	/*
	* Method: clear
	* Usage: map.clear();
	* -------------------
	* This method removes all entries from this map.
	*/
	void clear();

	/*
	* Method: put
	* Usage: map.put(key, value);
	* ---------------------------
	* This method associates key with value in this map. Any value
	* previously associated with this key is replaced by the new one.
	*/
	void put(int key, int value);

	/*
	* Method: get
	* Usage: value = map.get(key);
	* ----------------------------
	* If key is found in this map, this method returns the associated
	* value. If key is not found, the get mathod raises an error.
	* Clients can use the containsKey method to verify the presence
	* of a key in the map before attempting to get its value.
	*/
	bool get(int key, int *value);

	/*
	* Method: containsKey
	* Usage: if (map.containsKey(key))...
	* -----------------------------------
	* Returns true if there is an entry for key in this map,
	* false otherwise.
	*/
	bool containsKey(int key);

	/*
	* Method: remove
	* Usage: map.remove(key);
	* -----------------------
	* This method removes any entry for key from this map.
	* If there is no entry for the key, the map is unchanged.
	*/
	void remove(int key);

private:
	/* Constants */
	static const int INITIAL_SIZE = 101;

	/* Type for a linked list cell */
	typedef struct cellT 
	{
		int key;
		int value;
		cellT *link;
	} CellT;

	/* Instance variables */
	cellT **buckets; /* A dynamic array of the buckets */
	int nBuckets; /* Allocated size of the buckets array */
	int nEntries; /* The number of entries in the map */

	/* Private method prototypes */
	int hash(int key);
	cellT *findCell(cellT *chain, int key);
	void deleteChain(cellT *chain);
};
#endif
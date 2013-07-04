#ifndef __TRIPLE_H__
#define __TRIPLE_H__


//******************triple template*************************
//based on the STL design of the pair<> template
template <class T1, class T2, class T3>
struct triple{
	//type names for the values
	typedef T1 first_type;
	typedef T2 second_type;
	typedef T3 third_type;

	//member
	T1 first;
	T2 second;
	T3 third;

	/*default constructor
	*	T1() and T2() and T3() force initialization for built in types
	*/
	triple()
		:	first(T1()),second(T2()),third(T3())
	{
	}

	//constructor for 3 values
	triple(const T1 &a, const T2 &b, const T3 &c)
		:	first(a), second(b), third(c)
	{
	}

        // equality operator
        bool operator== (triple<T1,T2,T3> & t2)
        {
          return ((first == t2.first) && (second == t2.second) && (third == t2.third));
        }

        bool operator!= (triple<T1,T2,T3> & t2)
        {
          return ((first != t2.first) || (second != t2.second) || (third != t2.third));
        }

	//copy constructor with implicit conversions
	template<class U, class V, class W>
	triple(const triple<U,V,W> &t)
		:	first(t.first), second(t.second), third(t.third)
	{
	}

};


template <class T1, class T2, class T3>
class tripleCompare{
public:
  bool operator()(const triple<T1,T2,T3> & x, const triple<T1,T2,T3> & y) const
  {
    if (x.first < y.first)
      return true;

    if (y.first < x.first)
      return false;

    // now, we know that x.first == y.first
    // compare on second

    if (x.second < y.second)
      return true;

    if (y.second < x.second)
      return false;

    // now, we know that (x.first == y.first) && (x.second == y.second)
    // compare on second

    if (x.third < y.third)
      return true;

    return false;

  }
};

//**************quad template********************
//based on the design of the STL pair<> template
template <class T1, class T2, class T3, class T4>
struct quad{
	//type names for the values
	typedef T1 first_type;
	typedef T2 second_type;
	typedef T3 third_type;
	typedef T4 fourth_type;

	//member
	T1 first;
	T2 second;
	T3 third;
	T4 fourth;

	/*default constructor
	*	T1() and T2() and T3() force initialization for built in types
	*/
	quad()
		:	first(T1()),second(T2()),third(T3()),fourth(T4())
	{
	}

	//constructor for 3 values
	quad(const T1 &a, const T2 &b, const T3 &c, const T4 &d)
		:	first(a), second(b), third(c), fourth(d)
	{
	}

	//copy constructor with implicit conversions
	template<class U, class V, class W, class X>
	quad(const quad<U,V,W,X> &q)
		:	first(q.first), second(q.second), third(q.third), fourth(q.fourth)
	{
	}

};


#endif

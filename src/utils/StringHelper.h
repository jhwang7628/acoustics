#ifndef STRINGHELPER_H
#define STRINGHELPER_H


#include <vector>
#include <string> 
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/erase.hpp>


using namespace std; 


class StringHelper
{
    public: 
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 

        static vector<string> strtokall( const string & str, const char * delim )
        {
        
            boost::char_separator<char> sep(delim);
            tokenizer tokens(str, sep);

            vector<string> alltokens;
            copy(tokens.begin(), tokens.end(), back_inserter<vector<string> >(alltokens));
            //copy(tokens.begin(), tokens.end(), allfaces);
        
            return alltokens;
        }

        static vector<string> strtokall( const char * strchar, const char * delim )
        {

            string str( strchar );
            vector<string> alltokens = StringHelper::strtokall( string( str ), delim ); 
        
            //boost::char_separator<char> sep(delim);
            //tokenizer tokens(str, sep);

            //vector<string> alltokens;
            //copy(tokens.begin(), tokens.end(), back_inserter<vector<string> >(alltokens));
            //copy(tokens.begin(), tokens.end(), allfaces);
        
            return alltokens;
        }

        static bool matchSubstr( const string &str, const string &substr) 
        {
            boost::cmatch what; 
            boost::regex ex( substr.c_str() );
            return (boost::regex_match(str.c_str(), what, ex)); 
        }

        /* 
         * Wrapper for finding substring using boost regular expression 
         *
         * Example: 
         *   (impedance.txt)\.\?\$ search for impedance.txt at the end
         */
        static bool searchSubstr( const string & str, const string & substr )
        {

            boost::regex ex( substr.c_str() );
            boost::smatch what; 

            return (boost::regex_search( str, what, ex ) );
        }

        /* 
         * Wrapper for finding substrings using boost regular expression 
         */
        static bool searchSubstrs( const string & str, const vector<string> & substrs )
        {

            size_t count = 0; 
            for ( size_t ii=0; ii<substrs.size(); ii++ ) 
            {
                if (searchSubstr( str, substrs[ii] ))
                    count ++;

            }

            return (count == substrs.size() ? 1 : 0);
        }


        /* 
         * Wrapper for finding substrings using boost regular expression 
         */
        static bool searchSubstrs( const string & str, const string & substr1, const string & substr2 )
        {

            vector<string> tmp; 
            tmp.push_back( substr1 ); 
            tmp.push_back( substr2 ); 

            return searchSubstrs( str, tmp );


        }

        static string stripChar( const string & str, const char & chr )
        {
            string str_mod( str ); 
            
            boost::erase_all( str_mod, &chr ); 

            return str_mod; 
        }

};



#endif

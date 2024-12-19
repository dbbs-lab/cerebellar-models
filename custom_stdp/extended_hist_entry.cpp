#include "extended_hist_entry.h"
#include <iostream>
namespace nest
{

ExtendedHistEntry::ExtendedHistEntry( double t, double Kminus, double Kminus_triplet, double offset, size_t access_counter )
  : t_( t )
  , Kminus_( Kminus )
  , Kminus_triplet_( Kminus_triplet )
  , offset_ ( offset )
  , access_counter_( access_counter )
{
  std::cout << "Created extended history entry with offset = " << offset_ << std::endl;
}

}

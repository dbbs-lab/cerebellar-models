#ifndef EXTENDED_HIST_ENTRY
#define EXTENDED_HIST_ENTRY

#include <cmath>

namespace nest
{

class ExtendedHistEntry
{
public:
  ExtendedHistEntry( double t, double Kminus, double Kminus_triplet, double offset, size_t access_counter );

  double t_;              //!< point in time when spike occurred (in ms)
  double Kminus_;         //!< value of Kminus at that time
  double Kminus_triplet_; //!< value of triplet STDP Kminus at that time
  double offset_; //!< new field for sneakily storing the flag
  size_t access_counter_; //!< access counter to enable removal of the entry, once all neurons read it
};


}

#endif

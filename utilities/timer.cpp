#include "timer.h"
// #define DISABLE_TIMER

using namespace mytimes;

template< typename T >
std::string to_pretty_str( T n )
{
    std::stringstream ss;
    ss.precision(3);
    if( n<1e4 )
        ss << n;
    else if( n<1e6 )
    {
        int thousands = std::floor( n/1000 );
        ss << thousands << ","
           << std::setw(3) << std::setfill('0')
           << n-1000*thousands;
        }
    else if( n<1e9 )
        ss << n/1.e6 << " million";
    else if( n<1e12 )
        ss << n/1.e9 << " billion";
    else
        ss << n/1.e12 << " trillion";
    return ss.str();
}


timer::timer(const string &label) : label_(label) {}
            
timer::~timer() {
    #ifndef DISABLE_TIMER
        stop(); 
        if(!silent){
            print();
        }
    #endif
}

void timer::set_label(string l) {
    #ifndef DISABLE_TIMER
        label_ = l;
    #endif
}

const timer& mytimes::timer::print() const {
    #ifndef DISABLE_TIMER
        double count = total_time.count();
        cout << to_pretty_str(count)
             << " seconds spent in program unit '" << label_ << "'" 
             << endl;
    #endif

    return *this;
}

timer &mytimes::timer::start() {
    #ifndef DISABLE_TIMER
        start_time = system_clock::now();
    #endif

    return *this;
}

timer &timer::stop() {
    #ifndef DISABLE_TIMER
        end_time = system_clock::now();
        total_time += end_time - start_time;
    #endif

    return *this;
}
            
timer &timer::reset() {
    #ifndef DISABLE_TIMER
        start_time = system_clock::now();
        total_time = start_time - start_time;
    #endif

    return *this;
}

timer &timer::reset_and_print() {
    #ifndef DISABLE_TIMER
        print();
        reset();
    #endif

    return *this;
}

void timer::silence() {
    #ifndef DISABLE_TIMER
        silent = true;
    #endif
}

bool timer::silent = false;
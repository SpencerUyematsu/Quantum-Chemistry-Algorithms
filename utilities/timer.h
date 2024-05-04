#include <chrono>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>


namespace mytimes{
    using namespace std::chrono;
    using namespace std;

    class timer{
        private:
            string label_;
            system_clock::time_point start_time;
            system_clock::time_point end_time;
            duration<double> total_time;

        public:
            static bool silent;
            explicit timer(const string &label);
            ~timer();

            void set_label(string l);

            const timer &print() const;

            timer &start();
            timer &stop();
            timer &reset();

            timer &reset_and_print();

            static void silence();
    };
}